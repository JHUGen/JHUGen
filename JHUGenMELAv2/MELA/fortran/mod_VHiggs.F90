!--YaofuZhou-----------------------------------------
module ModVHiggs
  use ModParameters
  use ModMisc
  implicit none
  private


  !----- notation for subroutines
  public :: EvalAmp_VHiggs,EvalUnpolAmpSq_gg_VH,TreeAmp_qqb_ZH,TreeAmp_qqbg_ZH
  public :: Dipoles_qqb_VH,CalcFormFactor1
  real,parameter,public :: alpha_ii=1d0

contains

subroutine EvalAmp_VHiggs(id,helicity,MomExt,me2)
      integer, intent(in) :: id(9)
      real(8), intent(in) :: helicity(9)
      real(8), intent(in) :: MomExt(1:4,1:9)
      real(8), intent(out) :: me2
      real(8) :: mass(3:5,1:2)
      integer :: i
      complex(8) amplitude, A_VV(1:4), amptest
      integer :: idin(9)
      real(8) :: helin(9)
      real(8) :: pin(1:4,1:9)

      idin(:)=id(:)
      helin(:)=helicity(:)
      pin(:,:)=MomExt(:,:)
      if(id(2).eq.convertLHE(Pho_)) then
         call swap(idin(1),idin(2))
         call swap(helin(1),helin(2))
         call swap(pin(:,1),pin(:,2))
      endif
      if(id(7).eq.convertLHE(Pho_)) then
         call swap(idin(6),idin(7))
         call swap(helin(6),helin(7))
         call swap(pin(:,6),pin(:,7))
      endif

      do i=3,5
        mass(i,1) = getMass(convertLHEreverse(idin(i)))
        mass(i,2) = getDecayWidth(convertLHEreverse(idin(i)))
      enddo

      A_VV(:)=czero
      if(idin(1).ne.convertLHE(Pho_) .and. idin(6).ne.convertLHE(Pho_)) then
         !print *,"Case1:"
         A_VV(1)=MATRIXELEMENT0(pin,mass,helin,idin,(/.false., .false./))
         if(includeGammaStar) then
            !print *,"Case2:"
            A_VV(2) = MATRIXELEMENT0(pin,mass,helin,idin,(/.false., .true./))
            !print *,"Case3:"
            A_VV(3) = MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .false./))
            !print *,"Case4:"
            A_VV(4) = MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .true./))
         endif
      else if(idin(1).eq.convertLHE(Pho_) .and. idin(6).eq.convertLHE(Pho_)) then
         !print *,"Case5:"
         A_VV(1)=MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .true./))
      else if(idin(1).eq.convertLHE(Pho_)) then
         !print *,"Case6:"
         A_VV(1)=MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .false./))
         if(includeGammaStar) then
            !print *,"Case7:"
            A_VV(2) = MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .true./))
         endif
      else !if(idin(6).eq.convertLHE(Pho_)) then
         !print *,"Case8:"
         A_VV(1)=MATRIXELEMENT0(pin,mass,helin,idin,(/.false., .true./))
         if(includeGammaStar) then
            !print *,"Case9:"
            A_VV(2) = MATRIXELEMENT0(pin,mass,helin,idin,(/.true., .true./))
         endif
      endif
      amplitude = A_VV(1)+A_VV(2)+A_VV(3)+A_VV(4)

      ! XCHECK FROM DECAY ME
      !print *,pin(:,1)
      !print *,pin(:,2)
      !print *,pin(:,6)
      !print *,pin(:,7)
      !print *,idin(1)
      !print *,idin(2)
      !print *,idin(3)
      !print *,idin(4)
      !print *,idin(6)
      !print *,idin(7)
      !print *,helin(1)
      !print *,helin(2)
      !print *,helin(6)
      !print *,helin(7)
      !print *,amplitude
      !amptest = MATRIXELEMENT02(pin,mass,helin,idin)
      !print *,amptest
      !pause

      me2=dble(amplitude*dconjg(amplitude))
      return
end subroutine EvalAmp_VHiggs

!MATRIXELEMENT0.F
!VERSION 20160511
      complex(8) function MATRIXELEMENT0(MomExt,mass,helicity,id,useA)
      implicit none
      real(8), intent(in) :: MomExt(1:4,1:9)
      real(8), intent(in) :: mass(3:5,1:2)
      real(8), intent(in) :: helicity(9)
      integer, intent(in) :: id(9)
      logical, intent(in) :: useA(2)

      integer mu3,mu4
      real(8) qq,q3_q3,q4_q4,q5_q5
      complex(8) PVVX0P
      complex(8) Vcurrent1(4), Acurrent1(4), current1(4), currentVp1(4)
      complex(8) Vcurrent2(4), Acurrent2(4), current2(4), currentVp2(4)
      complex(8) Vpffcoupl(2,2)
      complex(8) POL1(3,4), POL2(3,4)
      complex(8) g_mu_nu(4,4), pp(4,4), epp(4,4)
      complex(8) PROP1, PROP2, PROP3
      complex(8) PROP_Vp1, PROP_Vp2
      complex(8) gFFZ, gFFA, gFFW
      complex(8) gVVP, gVVS1, gVVS2
      complex(8) ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      complex(8) gVVpP, gVVpS1, gVVpS2
      complex(8) ghzzp1_dyn,ghzzp2_dyn,ghzzp3_dyn,ghzzp4_dyn
      complex(8) gVpVP, gVpVS1, gVpVS2
      complex(8) ghzpz1_dyn,ghzpz2_dyn,ghzpz3_dyn,ghzpz4_dyn
      complex(8) gVpVpP, gVpVpS1, gVpVpS2
      complex(8) ghzpzp1_dyn,ghzpzp2_dyn,ghzpzp3_dyn,ghzpzp4_dyn

      MATRIXELEMENT0=czero
      if( &
         (id(1).ne.convertLHE(Pho_) .and. useA(1) .and. .not.includeGammaStar) .or. &
         (id(6).ne.convertLHE(Pho_) .and. useA(2) .and. .not.includeGammaStar) .or. &
         ((id(1)+id(2)).ne.0 .and. id(1).ne.convertLHE(Pho_) .and. useA(1)) .or. &
         ((id(6)+id(7)).ne.0 .and. id(6).ne.convertLHE(Pho_) .and. useA(2))      &
        ) then
         return
      endif

      Vcurrent1 = czero
      Acurrent1 = czero
      Vcurrent2 = czero
      Acurrent2 = czero
      Vpffcoupl=czero
      currentVp1=czero
      currentVp2=czero
      PROP_Vp1=czero
      PROP_Vp2=czero
      ghzzp1_dyn = czero
      ghzzp2_dyn = czero
      ghzzp3_dyn = czero
      ghzzp4_dyn = czero
      ghzpz1_dyn = czero
      ghzpz2_dyn = czero
      ghzpz3_dyn = czero
      ghzpz4_dyn = czero
      ghzpzp1_dyn = czero
      ghzpzp2_dyn = czero
      ghzpzp3_dyn = czero
      ghzpzp4_dyn = czero
      gVVpP = czero
      gVVpS1 = czero
      gVVpS2 = czero
      gVpVP = czero
      gVpVS1 = czero
      gVpVS2 = czero
      gVpVpP = czero
      gVpVpS1 = czero
      gVpVpS2 = czero

      gFFZ = ci*2d0*dsqrt(couplZffsq) ! = sqrt(gwsq/(1.0_dp-xw))
      gFFA = -ci*dsqrt(couplAffsq) ! = sqrt(gwsq*xw)
      gFFW = ci*dsqrt(couplWffsq) ! = sqrt(gwsq/2.0_dp)

      qq = -scr(MomExt(:,3),MomExt(:,4))
      if (id(1).eq.convertLHE(Pho_) .and. useA(1)) then
         q3_q3 = 0d0
      else
         q3_q3 = scr(MomExt(:,3),MomExt(:,3))
      endif
      if (id(6).eq.convertLHE(Pho_) .and. useA(2)) then
         q4_q4 = 0d0
      else
         q4_q4 = scr(MomExt(:,4),MomExt(:,4))
      endif
      q5_q5 = scr(MomExt(:,5),MomExt(:,5))
      PROP3 = PROPAGATOR(dsqrt(q5_q5),mass(5,1),mass(5,2))

      if (includeVprime) then
         if(.not.useA(1)) then
            Vpffcoupl(1,1)=GetVpffCoupling_VH(id(1), -1, ((id(1)+id(2)).ne.0))
            Vpffcoupl(1,2)=GetVpffCoupling_VH(id(1), +1, ((id(1)+id(2)).ne.0))
         endif
         if(.not.useA(2)) then
            Vpffcoupl(2,1)=GetVpffCoupling_VH(id(6), -1, ((id(6)+id(7)).ne.0))
            Vpffcoupl(2,2)=GetVpffCoupling_VH(id(6), +1, ((id(6)+id(7)).ne.0))
         endif
      endif


      if(.not.useA(1)) then
         PROP1 = PROPAGATOR(dsqrt(q3_q3),mass(3,1),mass(3,2))
         if(id(1).gt.0)then
            call FFV(id(2), MomExt(:,2), helicity(2), id(1), MomExt(:,1), helicity(1), Vcurrent1)
            call FFA(id(2), MomExt(:,2), helicity(2), id(1), MomExt(:,1), helicity(1), Acurrent1)
         else
            call FFV(id(1), MomExt(:,1), helicity(1), id(2), MomExt(:,2), helicity(2), Vcurrent1)
            call FFA(id(1), MomExt(:,1), helicity(1), id(2), MomExt(:,2), helicity(2), Acurrent1)
         endif

         ! Vpff current without the prefactor
         if (includeVprime) then
            if((id(1)*helicity(1)).le.0d0)then
               currentVp1=( &
                  Vcurrent1*(Vpffcoupl(1,1)+Vpffcoupl(1,2))*0.5 - &
                  Acurrent1*(Vpffcoupl(1,1)-Vpffcoupl(1,2))*0.5   &
                  )
            else
               currentVp1=( &
                  Vcurrent1*Vpffcoupl(1,2)                        &
                  )
            endif
         endif

         !WH
         if((id(1)+id(2)).ne.0)then
            if (includeVprime) then
               if (getMass(Wppr_).ge.0d0) then
                  !print *,"Compute prop for Wppr"
                  PROP_Vp1 = PROPAGATOR(dsqrt(q3_q3),getMass(Wppr_),getDecayWidth(Wppr_))
               else
                  PROP_Vp1 = PROPAGATOR(M_W,0d0,0d0)
               endif
               currentVp1 = currentVp1*gFFW*CKMbare(id(1),id(2))
            endif
            if((id(1)*helicity(1)).le.0d0)then
               current1=(Vcurrent1-Acurrent1)/2d0*gFFW*CKMbare(id(1),id(2))
            else
               current1=0d0
            endif
         !ZH
         else
            if (includeVprime) then
               if (getMass(Zpr_).ge.0d0) then
                  !print *,"Compute prop for Zpr"
                  PROP_Vp1 = PROPAGATOR(dsqrt(q3_q3),getMass(Zpr_),getDecayWidth(Zpr_))
               else
                  PROP_Vp1 = PROPAGATOR(M_Z,0d0,0d0)
               endif
               currentVp1 = currentVp1*gFFZ
            endif
            !e+ e- Z vertex for incoming states
            if((abs(id(1)).eq.11).or.(abs(id(1)).eq.13).or.(abs(id(1)).eq.15))then
              if((id(1)*helicity(1)).gt.0d0)then
                current1=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent1 -(0.5d0*T3lR)*Acurrent1
              else
                current1=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent1 -(0.5d0*T3lL)*Acurrent1
              endif
              current1=current1*gFFZ
            !u u~ Z vertex for incoming states
            else if((abs(id(1)).eq.2).or.(abs(id(1)).eq.4))then
              if((id(1)*helicity(1)).gt.0d0)then
                current1=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent1 -(0.5d0*T3uR)*Acurrent1
              else
                current1=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent1 -(0.5d0*T3uL)*Acurrent1
              endif
              current1=current1*gFFZ
            !d d~ Z vertex for incoming states
            else if((abs(id(1)).eq.1).or.(abs(id(1)).eq.3).or.(abs(id(1)).eq.5))then
              if((id(1)*helicity(1)).gt.0d0)then
                current1=(0.5d0*T3dR - QdR*sitW**2) *Vcurrent1 -(0.5d0*T3dR)*Acurrent1
              else
                current1=(0.5d0*T3dL - QdL*sitW**2) *Vcurrent1 -(0.5d0*T3dL)*Acurrent1
              endif
              current1=current1*gFFZ
            else
              current1=0d0
              currentVp1=0d0
              print *, "invalid incoming state"
            endif
         endif
      else
         if(abs(id(1)).eq.convertLHE(Pho_)) then
           PROP1=cone
           if((id(1)*helicity(1)).gt.0d0) then
             call POLARIZATION_SINGLE(MomExt(:,3),+1,Vcurrent1)
           else
             call POLARIZATION_SINGLE(MomExt(:,3),-1,Vcurrent1)
           endif
         else
           PROP1 = PROPAGATOR(dsqrt(q3_q3),0d0,0d0)
           if(id(1).gt.0)then
             call FFV(id(2), MomExt(:,2), helicity(2), id(1), MomExt(:,1), helicity(1), Vcurrent1)
           else
             call FFV(id(1), MomExt(:,1), helicity(1), id(2), MomExt(:,2), helicity(2), Vcurrent1)
           endif
         endif

         !ZH
         if(abs(id(1)).eq.convertLHE(Pho_))then
           current1=Vcurrent1
         !e+ e- Z vertex for incoming states
         else if((abs(id(1)).eq.11).or.(abs(id(1)).eq.13).or.(abs(id(1)).eq.15))then
           if((id(1)*helicity(1)).gt.0d0)then
             current1 = QlR*Vcurrent1
           else
             current1 = QlL*Vcurrent1
           endif
           current1=current1*gFFA
         !u u~ Z vertex for incoming states
         else if((abs(id(1)).eq.2).or.(abs(id(1)).eq.4))then
           if((id(1)*helicity(1)).gt.0d0)then
             current1 = QuR*Vcurrent1
           else
             current1 = QuL*Vcurrent1
           endif
           current1=current1*gFFA
         !d d~ Z vertex for incoming states
         else if((abs(id(1)).eq.1).or.(abs(id(1)).eq.3).or.(abs(id(1)).eq.5))then
           if((id(1)*helicity(1)).gt.0d0)then
             current1 = QdR*Vcurrent1
           else
             current1 = QdL*Vcurrent1
           endif
           current1=current1*gFFA
         else
           current1=0d0
           print *, "invalid incoming state"
         endif
      endif

      if(.not.useA(2)) then
         PROP2 = PROPAGATOR(dsqrt(q4_q4),mass(4,1),mass(4,2))

         if(id(6).gt.0)then
           call FFV(id(6), MomExt(:,6), helicity(6), id(7), MomExt(:,7), helicity(7), Vcurrent2)
           call FFA(id(6), MomExt(:,6), helicity(6), id(7), MomExt(:,7), helicity(7), Acurrent2)
         else
           call FFV(id(7), MomExt(:,7), helicity(7), id(6), MomExt(:,6), helicity(6), Vcurrent2)
           call FFA(id(7), MomExt(:,7), helicity(7), id(6), MomExt(:,6), helicity(6), Acurrent2)
         endif

         ! Vpff current without the prefactor
         if (includeVprime) then
            if((id(6)*helicity(6)).le.0d0)then
               currentVp2=( &
                  Vcurrent2*(Vpffcoupl(2,1)+Vpffcoupl(2,2))*0.5 - &
                  Acurrent2*(Vpffcoupl(2,1)-Vpffcoupl(2,2))*0.5   &
                  )
            else
               currentVp2=( &
                  Vcurrent2*Vpffcoupl(2,2)                        &
                  )
            endif
         endif

         !WH
         if((id(6)+id(7)).ne.0)then
           if (includeVprime) then
             if (getMass(Wppr_).ge.0d0) then
               PROP_Vp2 = PROPAGATOR(dsqrt(q4_q4),getMass(Wppr_),getDecayWidth(Wppr_))
             else
               PROP_Vp2 = PROPAGATOR(M_W,0d0,0d0)
             endif
             currentVp2 = currentVp2*gFFW*CKMbare(id(6),id(7))
           endif
           if((id(6)*helicity(6)).le.0d0)then
             current2=(Vcurrent2-Acurrent2)/2d0*gFFW*CKM(id(6),id(7))
           else
             current2=0d0
           endif
         !ZH
         else
           if (includeVprime) then
             if (getMass(Zpr_).ge.0d0) then
               PROP_Vp2 = PROPAGATOR(dsqrt(q4_q4),getMass(Zpr_),getDecayWidth(Zpr_))
             else
               PROP_Vp2 = PROPAGATOR(M_Z,0d0,0d0)
             endif
             currentVp2 = currentVp2*gFFZ
           endif
           !l+ l- Z vertex for final state
           if((abs(id(6)).eq.11).or.(abs(id(6)).eq.13))then
             if((id(6)*helicity(6)).gt.0d0)then
               current2=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent2 -(0.5d0*T3lR)*Acurrent2
             else
               current2=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent2 -(0.5d0*T3lL)*Acurrent2
             endif
             current2=current2*gFFZ*dsqrt(scale_alpha_Z_ll)
           !tau+ tau- Z vertex for final state
           else if((abs(id(6)).eq.15))then
             if((id(6)*helicity(6)).gt.0d0)then
               current2=(0.5d0*T3lR - QlR*sitW**2) *Vcurrent2 -(0.5d0*T3lR)*Acurrent2
             else
               current2=(0.5d0*T3lL - QlL*sitW**2) *Vcurrent2 -(0.5d0*T3lL)*Acurrent2
             endif
             current2=current2*gFFZ*dsqrt(scale_alpha_Z_tt)
           !u u~ Z vertex for final state
           else if((abs(id(6)).eq.2).or.(abs(id(6)).eq.4))then
             if((id(6)*helicity(6)).gt.0d0)then
               current2=(0.5d0*T3uR - QuR*sitW**2) *Vcurrent2 -(0.5d0*T3uR)*Acurrent2
             else
               current2=(0.5d0*T3uL - QuL*sitW**2) *Vcurrent2 -(0.5d0*T3uL)*Acurrent2
             endif
             current2=current2*gFFZ*dsqrt(scale_alpha_Z_uu)
           !d d~ Z vertex for final state
           else if((abs(id(6)).eq.1).or.(abs(id(6)).eq.3).or.(abs(id(6)).eq.5))then
             if((id(6)*helicity(6)).gt.0d0)then
               current2=(0.5d0*T3dR - QdR*sitW**2) *Vcurrent2 -(0.5d0*T3dR)*Acurrent2
             else
               current2=(0.5d0*T3dL - QdL*sitW**2) *Vcurrent2 -(0.5d0*T3dL)*Acurrent2
             endif
             current2=current2*gFFZ*dsqrt(scale_alpha_Z_dd)
           !nu nu~ Z vertex for final state
           else if((abs(id(6)).eq.12).or.(abs(id(6)).eq.14).or.(abs(id(6)).eq.16))then
             current2=(0.5d0*T3nL - QnL*sitW**2) *Vcurrent2 -(0.5d0*T3nL)*Acurrent2
             current2=current2*gFFZ*dsqrt(scale_alpha_Z_nn)
           else
             current2=0d0
             currentVp2 = 0d0
             print *, "invalid final state", id(6:7), helicity(6:7)
             stop
           endif
         endif
      else
         if(abs(id(6)).eq.convertLHE(Pho_)) then
           PROP2=cone
           if((id(6)*helicity(6)).gt.0d0) then
             call POLARIZATION_SINGLE(MomExt(:,4),+1,Vcurrent2)
           else
             call POLARIZATION_SINGLE(MomExt(:,4),-1,Vcurrent2)
           endif
           Vcurrent2 = dconjg(Vcurrent2)
         else
           PROP2 = PROPAGATOR(dsqrt(q4_q4),0d0,0d0)
           if(id(6).gt.0)then
             call FFV(id(6), MomExt(:,6), helicity(6), id(7), MomExt(:,7), helicity(7), Vcurrent2)
           else
             call FFV(id(7), MomExt(:,7), helicity(7), id(6), MomExt(:,6), helicity(6), Vcurrent2)
           endif
         endif

         !ZH
         if(abs(id(6)).eq.convertLHE(Pho_)) then
           current2=Vcurrent2
         !l+ l- Z vertex for final state
         else if((abs(id(6)).eq.11).or.(abs(id(6)).eq.13))then
           if((id(6)*helicity(6)).gt.0d0)then
             current2=QlR*Vcurrent2
           else
             current2=QlL*Vcurrent2
           endif
           current2=current2*gFFA*dsqrt(scale_alpha_Z_ll)
         !tau+ tau- Z vertex for final state
         else if((abs(id(6)).eq.15))then
           if((id(6)*helicity(6)).gt.0d0)then
             current2=QlR*Vcurrent2
           else
             current2=QlL*Vcurrent2
           endif
           current2=current2*gFFA*dsqrt(scale_alpha_Z_tt)
         !u u~ Z vertex for final state
         else if((abs(id(6)).eq.2).or.(abs(id(6)).eq.4))then
           if((id(6)*helicity(6)).gt.0d0)then
             current2=QuR*Vcurrent2
           else
             current2=QuL*Vcurrent2
           endif
           current2=current2*gFFA*dsqrt(scale_alpha_Z_uu)
         !d d~ Z vertex for final state
         else if((abs(id(6)).eq.1).or.(abs(id(6)).eq.3).or.(abs(id(6)).eq.5))then
           if((id(6)*helicity(6)).gt.0d0)then
             current2=QdR*Vcurrent2
           else
             current2=QdL*Vcurrent2
           endif
           current2=current2*gFFA*dsqrt(scale_alpha_Z_dd)
         !nu nu~ Z vertex for final state
         else if((abs(id(6)).eq.12).or.(abs(id(6)).eq.14).or.(abs(id(6)).eq.16))then
           current2=QnL*Vcurrent2
           current2=current2*gFFA*dsqrt(scale_alpha_Z_nn)
         else
           current2=0d0
           print *, "invalid final state", id(6:7), helicity(6:7)
           stop
         endif
      endif

      if(.not.(useA(1) .and. abs(id(1)).eq.convertLHE(Pho_))) then
         current1 = -current1 + scrc(MomExt(:,3),current1)/q3_q3
         currentVp1 = -currentVp1 + scrc(MomExt(:,3),currentVp1)/q3_q3
      endif
      if(.not.(useA(2) .and. abs(id(6)).eq.convertLHE(Pho_))) then
         current2 = -current2 + scrc(MomExt(:,4),current2)/q4_q4
         currentVp2 = -currentVp2 + scrc(MomExt(:,4),currentVp2)/q4_q4
      endif

      !print *,"current1=",current1
      !print *,"currentVp1=",currentVp1
      !print *,"PROP1=",PROP1
      !print *,"PROP_Vp1=",PROP_Vp1
      !print *,"current2=",current2
      !print *,"currentVp2=",currentVp2
      !print *,"PROP2=",PROP2
      !print *,"PROP_Vp2=",PROP_Vp2

      current1 = current1*PROP1
      current2 = current2*PROP2
      currentVp1 = currentVp1*PROP_Vp1
      currentVp2 = currentVp2*PROP_Vp2

!XVV vertex
      if(id(3).eq.convertLHE(Wp_))then
         call swap(q3_q3,q4_q4)
         call swap(current1,current2)
         call swap(currentVp1,currentVp2)
      endif

      if(.not.useA(1) .and. .not.useA(2)) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(1,q3_q3,q4_q4,q5_q5)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(2,q3_q3,q4_q4,q5_q5)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(3,q3_q3,q4_q4,q5_q5)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(4,q3_q3,q4_q4,q5_q5)

         if (includeVprime) then
            ghzzp1_dyn = HVVSpinZeroDynamicCoupling(12,q3_q3,q4_q4,q5_q5)
            ghzzp2_dyn = HVVSpinZeroDynamicCoupling(13,q3_q3,q4_q4,q5_q5)
            ghzzp3_dyn = HVVSpinZeroDynamicCoupling(14,q3_q3,q4_q4,q5_q5)
            ghzzp4_dyn = HVVSpinZeroDynamicCoupling(15,q3_q3,q4_q4,q5_q5)

            ghzpz1_dyn = HVVSpinZeroDynamicCoupling(12,q4_q4,q3_q3,q5_q5)
            ghzpz2_dyn = HVVSpinZeroDynamicCoupling(13,q4_q4,q3_q3,q5_q5)
            ghzpz3_dyn = HVVSpinZeroDynamicCoupling(14,q4_q4,q3_q3,q5_q5)
            ghzpz4_dyn = HVVSpinZeroDynamicCoupling(15,q4_q4,q3_q3,q5_q5)

            ghzpzp1_dyn = HVVSpinZeroDynamicCoupling(16,q3_q3,q4_q4,q5_q5)
            ghzpzp2_dyn = HVVSpinZeroDynamicCoupling(17,q3_q3,q4_q4,q5_q5)
            ghzpzp3_dyn = HVVSpinZeroDynamicCoupling(18,q3_q3,q4_q4,q5_q5)
            ghzpzp4_dyn = HVVSpinZeroDynamicCoupling(19,q3_q3,q4_q4,q5_q5)
         endif
      else if(useA(1) .and. useA(2)) then
         ghz1_dyn = czero
         ghz2_dyn = HVVSpinZeroDynamicCoupling(9,q3_q3,q4_q4,q5_q5)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(10,q3_q3,q4_q4,q5_q5)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(11,q3_q3,q4_q4,q5_q5)
      else if(useA(1)) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q3_q3,q5_q5)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q3_q3,q5_q5)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q3_q3,q5_q5)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q3_q3,q5_q5)
      else !if(useA(2)) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q5_q5)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q5_q5)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q5_q5)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q5_q5)
      endif

      gVVS1 = ghz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda**2 )
      gVVS2 = -( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda**2 )
      gVVP = -2d0*ghz4_dyn

      if(.not.useA(1) .and. .not.useA(2)) then
         if (includeVprime) then
            gVVpS1 = ghzzp1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpS2 = -( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpP = -2d0*ghzzp4_dyn

            gVpVS1 = ghzpz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVS2 = -( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVP = -2d0*ghzpz4_dyn

            gVpVpS1 = ghzpzp1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzpzp2_dyn + ghzpzp3_dyn*qq/Lambda**2 )
            gVpVpS2 = -( 2d0*ghzpzp2_dyn + ghzpzp3_dyn*qq/Lambda**2 )
            gVpVpP = -2d0*ghzpzp4_dyn
         endif
      endif


      call VVS1(g_mu_nu)
      call VVS2(MomExt(:,5),MomExt(:,5),pp)
      if(id(3).eq.convertLHE(Wp_))then
         call VVP(MomExt(:,4),-MomExt(:,3),epp)
      else
         call VVP(-MomExt(:,3),MomExt(:,4),epp)
      endif

! assemble everything and get iM
      MATRIXELEMENT0=(0d0,0d0)
      do mu3=1,4
      do mu4=1,4
         MATRIXELEMENT0 = MATRIXELEMENT0 +      &
            current1(mu3)*current2(mu4)*(       &
               gVVS1*g_mu_nu(mu3,mu4)         + &
               gVVS2*pp     (mu3,mu4)         + &
               gVVP *epp    (mu3,mu4)           &
            )                                   &
            +                                   &
            current1(mu3)*currentVp2(mu4)*(    &
               gVVpS1*g_mu_nu(mu3,mu4)        + &
               gVVpS2*pp     (mu3,mu4)        + &
               gVVpP *epp    (mu3,mu4)          &
            )                                   &
            +                                   &
            currentVp1(mu3)*current2(mu4)*(    &
               gVpVS1*g_mu_nu(mu3,mu4)        + &
               gVpVS2*pp     (mu3,mu4)        + &
               gVpVP *epp    (mu3,mu4)          &
            )                                   &
            +                                   &
            currentVp1(mu3)*currentVp2(mu4)*( &
               gVpVpS1*g_mu_nu(mu3,mu4)       + &
               gVpVpS2*pp     (mu3,mu4)       + &
               gVpVpP *epp    (mu3,mu4)         &
            )
      enddo !mu4
      enddo !mu3
      MATRIXELEMENT0 = MATRIXELEMENT0*ci/vev

      if(H_DK.eqv..false.)then
        MATRIXELEMENT0=MATRIXELEMENT0 *PROP3
      else if(id(8).ne.Not_a_particle_) then
        MATRIXELEMENT0=MATRIXELEMENT0 *PROP3 &
        *(kappa*FFS(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)) &
         +kappa_tilde*FFP(id(8), MomExt(:,8), helicity(8), id(9), MomExt(:,9), helicity(9)))&
        *(-ci/vev*getMass(convertLHEreverse(id(8))))
      else
        MATRIXELEMENT0=czero
      endif

      !print *,"MATRIXELEMENT0=",MATRIXELEMENT0

      return
      END function



function GetVpffCoupling_VH(pdgid, hel, useWp)
integer, intent(in) :: pdgid
integer, intent(in) :: hel
logical, intent(in) :: useWp
complex(8) :: GetVpffCoupling_VH
   GetVpffCoupling_VH=VpffCoupling_PDG(pdgid,hel,useWp)
   if(useWp) then
     GetVpffCoupling_VH = GetVpffCoupling_VH / bL ! Bc of the couplings formalism from decay
   else
     GetVpffCoupling_VH = GetVpffCoupling_VH / 2.0_dp
   endif
end function


!ANGLES.F
!VERSION 20130531
!
!A subroutine that calculates the polar and azimuthal angles of a
!given vector(4) in terms of their sin and cos, which will be
!returned by the array sincos(4).
      subroutine ANGLES(sincos, vector)
      implicit none
!     real(8) Pi
      real(8) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      real(8) sincos(4), vector(4), abs3p, phi
!sincos(1)=cos(theta)
!sincos(2)=sin(theta)
!sincos(3)=cos(phi)
!sincos(4)=sin(phi)

!|3-momentum|
      abs3p = dsqrt(vector(2)**2+vector(3)**2+vector(4)**2)

!if |3-momentum|=0
      if(abs3p.lt.epsilon)then
        sincos(1)=1d0
        sincos(2)=0d0
      else
        sincos(1)=vector(4)/abs3p
        sincos(2)=dsqrt((1d0+sincos(1))*(1d0-sincos(1)))
      endif

!if colinear
      if(dabs(vector(3)).lt.epsilon)then
        phi=0d0
      else
        if(dabs(vector(2)).lt.epsilon)then
           phi=(TwoPi/2d0)/2d0 * dsign(1d0,vector(3))
        else
           phi=datan(vector(3)/vector(2))
        endif
      endif
!shift phi so that 0 < phi < 2Pi
      if(vector(2).lt.0d0)then
         phi=phi+Pi
      endif
      if(phi.lt.0d0)then
         phi=phi+Twopi
      endif
!     print *,phi
      sincos(3)=dcos(phi)
      sincos(4)=dsin(phi)

      return
      END subroutine ANGLES

!ANTISYMMETRIC2.F
!VERSION 20130702
      subroutine ANTISYMMETRIC2(p1,p2,epp)

      implicit none
!     real(8) Pi
      real(8) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision

      complex(8) p1(4), p2(4)
      complex(8) epp(4,4)
!      real(8) ANTISYMMETRIC
      integer i,j,k,l

!      external ANTISYMMETRIC

!     do i=1,4
!     do j=1,4
      epp(i,j)=0d0
!     enddo
!     enddo

      do i=1,4
      do j=1,4
      do k=1,4
      do l=1,4
          epp(i,j)=epp(i,j)+ANTISYMMETRIC(i,j,k,l)*p1(k)*p2(l)
      enddo
      enddo
      enddo
      enddo

      return
      END subroutine ANTISYMMETRIC2

!ANTISYMMETRIC.F
!VERSION 20130618

!returns the element of the rank-4 COVARIANT total antysymmetric
!tensor.
!ANTISYMMETRIC(0,1,2,3)=1.
      real(8) function ANTISYMMETRIC(i,j,k,l)

      implicit none
!     include '../COMMON.INI'

      integer i,j,k,l

      ANTISYMMETRIC=dble((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l))/12d0

      return
      END function ANTISYMMETRIC

!CONTRA_FIELD_TENSOR.F
!VERSION 20130630

!in T_mu_nu(lambda, ^mu, ^nu) returns the contrvariant field tensors
!for given polarization vectors.
      subroutine CONTRA_FIELD_TENSOR(POL, T_mu_nu)

      implicit none
!     real(8) Pi
      real(8) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      complex(8) epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      complex(8) epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      complex(8) POL(3,4), T_mu_nu(5,4,4)

      call CONTRA_OUTER(POL(1,:), POL(1,:), epep)
      call CONTRA_OUTER(POL(2,:), POL(2,:), emem)
      call CONTRA_OUTER(POL(1,:), POL(3,:), epe0)
      call CONTRA_OUTER(POL(3,:), POL(1,:), e0ep)
      call CONTRA_OUTER(POL(2,:), POL(3,:), eme0)
      call CONTRA_OUTER(POL(3,:), POL(2,:), e0em)
      call CONTRA_OUTER(POL(3,:), POL(3,:), e0e0)
      call CONTRA_OUTER(POL(1,:), POL(2,:), epem)
      call CONTRA_OUTER(POL(2,:), POL(1,:), emep)

!lambda = +2
      T_mu_nu(1,:,:)=epep
!lambda = -2
      T_mu_nu(2,:,:)=emem
!lambda = +3
      T_mu_nu(3,:,:)=(epe0+e0ep)/dsqrt(2d0)
!lambda = -3
      T_mu_nu(4,:,:)=(eme0+e0em)/dsqrt(2d0)
!lambda = 0
      T_mu_nu(5,:,:)=(epem+emep)/dsqrt(6d0) + e0e0/dsqrt(1.5d0)


      return
      END subroutine CONTRA_FIELD_TENSOR

!CONTRA_OUTER.F
!VERSION 20130620

!in pp(_mu, _nu) returns the CONTRAVARIANT tensor of
!p1 p2 outer product.
      subroutine CONTRA_OUTER(p1,p2,pp)

      implicit none
!     include '../COMMON.INI'
      complex(8) p1(4), p2(4)
      complex(8) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
        enddo
      enddo

      return
      END subroutine CONTRA_OUTER

!COVARIANT_FIELD_TENSOR.F
!VERSION 20130630

!in T_mu_nu(lambda, _mu, _nu) returns the covariant field tensors
!for given polarization vectors.
      subroutine COVARIANT_FIELD_TENSOR(POL, T_mu_nu)

      implicit none
!     real(8) Pi
      real(8) Twopi, Fourpi, epsilon
!     parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twopi = 2d0 * Pi )
      parameter( Fourpi = 4d0 * Pi )
      parameter( epsilon = 1d-13 ) !a small quantity slightly above machine precision
      complex(8) epep(4,4),emem(4,4),epe0(4,4),eme0(4,4),e0e0(4,4)
      complex(8) epem(4,4),e0ep(4,4),e0em(4,4),emep(4,4)
      complex(8) POL(3,4), T_mu_nu(5,4,4)

      call COVARIANT_OUTER(POL(1,:), POL(1,:), epep)
      call COVARIANT_OUTER(POL(2,:), POL(2,:), emem)
      call COVARIANT_OUTER(POL(1,:), POL(3,:), epe0)
      call COVARIANT_OUTER(POL(3,:), POL(1,:), e0ep)
      call COVARIANT_OUTER(POL(2,:), POL(3,:), eme0)
      call COVARIANT_OUTER(POL(3,:), POL(2,:), e0em)
      call COVARIANT_OUTER(POL(3,:), POL(3,:), e0e0)
      call COVARIANT_OUTER(POL(1,:), POL(2,:), epem)
      call COVARIANT_OUTER(POL(2,:), POL(1,:), emep)

!lambda = +2
      T_mu_nu(1,:,:)=epep
!lambda = -2
      T_mu_nu(2,:,:)=emem
!lambda = +3
      T_mu_nu(3,:,:)=(epe0+e0ep)/dsqrt(2d0)
!lambda = -3
      T_mu_nu(4,:,:)=(eme0+e0em)/dsqrt(2d0)
!lambda = 0
      T_mu_nu(5,:,:)=(epem+emep)/dsqrt(6d0) + e0e0/dsqrt(1.5d0)

      return
      END subroutine COVARIANT_FIELD_TENSOR

!COVARIANT_OUTER.F
!VERSION 20130620

!in pp(_mu, _nu) returns the COVARIANT tensor of p1 p2
!outer product.
      subroutine COVARIANT_OUTER(p1,p2,pp)

      implicit none
!     include '../COMMON.INI'
      complex(8) p1(4), p2(4)
      complex(8) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
          if( ( (mu.ne.1).and.(nu.eq.1) ).or. &
              ( (mu.eq.1).and.(nu.ne.1) ) )then
            pp(mu,nu)=-pp(mu,nu)
          endif
        enddo
      enddo

      return
      END subroutine COVARIANT_OUTER

!COVARIANT_VECTOR.F
!VERSION 20130703

!returns the component of the COVARIANT vector for given 4-vector
!and Lorentz index.
      complex(8) function COVARIANT_VECTOR(p,mu)

      implicit none

      complex(8) p(4)
      integer mu

      if(mu.ne.1)then
        COVARIANT_VECTOR = -p(mu)
      else
        COVARIANT_VECTOR = p(mu)
      endif

      return
      END function COVARIANT_VECTOR

!FFP.A
!VERSION 20130522

!returns i.Psi~(p1,s1).gamma5.Psi(p2,s2) for massless states
      complex(8) function FFP(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

      real(8) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(8) sqrt_pp1Dpp2

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
        FFP=0d0

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
        FFP=0d0
      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
        FFP=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
        FFP= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
      else
        sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
        FFP=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
      endif

      FFP=FFP*(0d0,-1d0)

      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        FFP=-dconjg(FFP)

      endif

      return
      END function FFP

!FFA.F
!VERSION 20130523

!in Acurrent(4) returns Psi~(p1,s1).gamma^mu.gamma5.Psi(p2,s2) for massless
!states.
      subroutine FFA(pdg_code1, p1, h1, pdg_code2, p2, h2, Acurrent)

      implicit none
      real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

      real(8) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(8) sqrt_pp1Dpp2, sqrt_pp1Xpp2
      complex(8) Acurrent(4)
      integer mu

      Acurrent = (0d0,0d0)

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).lt.0d0)then
        do mu=1,4
          Acurrent(mu)=0d0
        enddo

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. &
               ( dabs( p2(1)+p2(4) ).lt.epsilon ) )then

        Acurrent(1)= 2d0*dsqrt(p1(1)*p2(1))
        Acurrent(2)= 0d0
        Acurrent(3)= 0d0
        Acurrent(4)=-Acurrent(1)

      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then

        Acurrent(1)= dsqrt( 2d0*p1(1) / ( p2(1)+p2(4) ) ) &
                    *( p2(2) + (0d0,1d0)*p2(3) )
        Acurrent(2)= dsqrt( 2d0*p1(1) * ( p2(1)+p2(4) ) )
        Acurrent(3)= (0d0,1d0)*Acurrent(2)
        Acurrent(4)=-Acurrent(1)

      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then

        Acurrent(1)= dsqrt( 2d0*p2(1) / ( p1(1)+p1(4) ) ) &
                    *( p1(2) - (0d0,1d0)*p1(3) )
        Acurrent(2)= dsqrt( 2d0*p2(1) * ( p1(1)+p1(4) ) )
        Acurrent(3)=-(0d0,1d0)*Acurrent(2)
        Acurrent(4)=-Acurrent(1)

      else

        sqrt_pp1Dpp2= dsqrt( (p1(1)+p1(4)) / (p2(1)+p2(4)) )
        sqrt_pp1Xpp2= dsqrt( (p1(1)+p1(4)) * (p2(1)+p2(4)) )
        Acurrent(1)= sqrt_pp1Xpp2         &
                    +( p1(2) - (0d0,1d0)*p1(3) )     &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
        Acurrent(2)= ( p1(2) - (0d0,1d0)*p1(3) )/sqrt_pp1Dpp2   &
                    +( p2(2) + (0d0,1d0)*p2(3) )*sqrt_pp1Dpp2
        Acurrent(3)= ( (0d0,1d0)*p1(2) + p1(3) )/sqrt_pp1Dpp2   &
                    -( (0d0,1d0)*p2(2) - p2(3) )*sqrt_pp1Dpp2
        Acurrent(4)=sqrt_pp1Xpp2          &
                    -( p1(2) - (0d0,1d0)*p1(3) )                &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
      endif

!     print *, pdg_code1,h1,dble(pdg_code1)*h1,'!'
      if( (dble(pdg_code1)*h1) .lt. 0d0)then
!       print *, Acurrent
        do mu=1,4
          Acurrent(mu)=-dconjg(Acurrent(mu))
        enddo
      endif

      return
      END subroutine FFA

!FFS.F
!VERSION 20130522

!returns Psi~(p1,s1).Psi(p2,s2) for massless states
      complex(8) function FFS(pdg_code1, p1, h1, pdg_code2, p2, h2)

      implicit none
      real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

      real(8) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(8) sqrt_pp1Dpp2

      FFS = (0d0,0d0)

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).gt.0d0)then
        FFS=0d0

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. (   dabs( p2(1)+p2(4) ).lt.epsilon ) )then
        FFS=0d0
      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then
        FFS=-dsqrt(2d0*p1(1)*(p2(1)+p2(4)))
      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then
        FFS= dsqrt(2d0*p2(1)*(p1(1)+p1(4)))
      else
        sqrt_pp1Dpp2 = dsqrt((p1(1)+p1(4))/(p2(1)+p2(4)))
        FFS=(p2(2)-(0d0,1d0)*p2(3))*sqrt_pp1Dpp2- (p1(2)-(0d0,1d0)*p1(3))/sqrt_pp1Dpp2
      endif

      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        FFS=-dconjg(FFS)

      endif

      return
      END function FFS

!FFV.F
!VERSION 20130523

!in Vcurrent(4) returns Psi~(p1,s1).gamma^mu.Psi(p2,s2) for massless
!states.
      subroutine FFV(pdg_code1, p1, h1, pdg_code2, p2, h2, Vcurrent)

      implicit none
      real(8), parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

      real(8) p1(4), p2(4), h1, h2
      integer pdg_code1, pdg_code2
      real(8) sqrt_pp1Dpp2, sqrt_pp1Xpp2
      complex(8) Vcurrent(4)
      integer mu

      Vcurrent = (0d0,0d0)

      if( ( dble(pdg_code1) *h1* dble(pdg_code2) *h2 ).lt.0d0)then
        do mu=1,4
          Vcurrent(mu)=0d0
        enddo

      else if( ( dabs( p1(1)+p1(4) ).lt.epsilon ).and. ( dabs( p2(1)+p2(4) ).lt.epsilon ) )then

        Vcurrent(1)= 2d0*dsqrt(p1(1)*p2(1))
        Vcurrent(2)= 0d0
        Vcurrent(3)= 0d0
        Vcurrent(4)=-Vcurrent(1)

      else if(   dabs( p1(1)+p1(4) ).lt.epsilon )then

        Vcurrent(1)= dsqrt( 2d0*p1(1) / ( p2(1)+p2(4) ) ) *( p2(2) + (0d0,1d0)*p2(3) )
        Vcurrent(2)= dsqrt( 2d0*p1(1) * ( p2(1)+p2(4) ) )
        Vcurrent(3)= (0d0,1d0)*Vcurrent(2)
        Vcurrent(4)=-Vcurrent(1)

      else if(   dabs( p2(1)+p2(4) ).lt.epsilon )then

        Vcurrent(1)= dsqrt( 2d0*p2(1) / ( p1(1)+p1(4) ) ) *( p1(2) - (0d0,1d0)*p1(3) )
        Vcurrent(2)= dsqrt( 2d0*p2(1) * ( p1(1)+p1(4) ) )
        Vcurrent(3)=-(0d0,1d0)*Vcurrent(2)
        Vcurrent(4)=-Vcurrent(1)

      else

        sqrt_pp1Dpp2= dsqrt( (p1(1)+p1(4)) / (p2(1)+p2(4)) )
        sqrt_pp1Xpp2= dsqrt( (p1(1)+p1(4)) * (p2(1)+p2(4)) )
        Vcurrent(1)= sqrt_pp1Xpp2  &
                    +( p1(2) - (0d0,1d0)*p1(3) )   &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
        Vcurrent(2)= ( p1(2) - (0d0,1d0)*p1(3) )/sqrt_pp1Dpp2  &
                    +( p2(2) + (0d0,1d0)*p2(3) )*sqrt_pp1Dpp2
        Vcurrent(3)= ( (0d0,1d0)*p1(2) + p1(3) )/sqrt_pp1Dpp2  &
                    -( (0d0,1d0)*p2(2) - p2(3) )*sqrt_pp1Dpp2
        Vcurrent(4)=sqrt_pp1Xpp2  &
                    -( p1(2) - (0d0,1d0)*p1(3) )   &
                    *( p2(2) + (0d0,1d0)*p2(3) )/sqrt_pp1Xpp2
      endif

      if( (dble(pdg_code1)*h1) .lt. 0d0)then
        do mu=1,4
          Vcurrent(mu)=dconjg(Vcurrent(mu))
        enddo
      endif

      return
      END subroutine FFV


!INV_LORENTZ.F
!VERSION 20130602
!
!A subroutine that performs a general inverse boost to a four vector
!(vector) based on another four vector (boost). The primed and
!unprimed frames have their axes in parallel to one another.
!Rotation is not performed by this subroutine.
      subroutine INV_LORENTZ(vector, boost)

      implicit none

      real(8) vector(4), boost(4)
      real(8) lambda(4,4), vector_copy(4)
      real(8) beta(2:4), beta_sq, gamma
      integer i,j

      do i=2,4
        beta(i) = -boost(i)/boost(1)
      enddo

      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2

      gamma = 1d0/dsqrt(1d0-beta_sq)

      lambda(1,1) = gamma

      do i=2,4
        lambda(1,i) = gamma*beta(i)
        lambda(i,1) = lambda(1,i)
      enddo

      do i=2,4
      do j=2,4
        lambda(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
      enddo
      enddo

!apply boost to vector1
      vector_copy = vector
      vector = 0d0
      do i=1,4
      do j=1,4
        vector(i) = vector(i) + lambda(i,j)*vector_copy(j)
      enddo
      enddo

      return
      END subroutine INV_LORENTZ

!KRONECKER_DELTA.F

!KRONECKER_DELTA(i,j)
!A function that returns 1 if i=j, and 0 otherwise.
      real(8) function KRONECKER_DELTA(i,j)
      integer i,j
      if(i.eq.j)then
        KRONECKER_DELTA = 1d0
      else
        KRONECKER_DELTA = 0d0
      endif

      return
      end function KRONECKER_DELTA

!METRIC.F
!VERSION 20130524

!in METRIC returns the element of the Minkovski metric, with
!signature (1,-1,-1,-1), for given (_mu, _nu) or given (^mu, ^nu).
      real(8) function METRIC(mu,nu)

      implicit none
      integer mu, nu

      if(mu.ne.nu)then
        METRIC=0d0
      else if(mu.eq.1)then
        METRIC=1d0
      else
        METRIC=-1d0
      endif

      return
      END function METRIC

!`.F
!VERSION 20130524

!in POL(lambda,^mu) returns the polarization vectors for given
!4-momentum p
      subroutine POLARIZATION(p, POL)

      implicit none
      real(8) p(4)
      complex(8) POL(3,4)

      call POLARIZATION_SINGLE(p,+1,POL(1,:))
      call POLARIZATION_SINGLE(p,-1,POL(2,:))
      call POLARIZATION_SINGLE(p, 0,POL(3,:))

      return
      END subroutine POLARIZATION

      subroutine POLARIZATION_SINGLE(p, lambda, POL)

      implicit none
      real(8) p(4), sincos(4), inv_mass, abs3p
      complex(8) POL(4)
      integer lambda

      POL(:)=czero

      call ANGLES(sincos, p)
      !sincos(1)=cos(theta)
      !sincos(2)=sin(theta)
      !sincos(3)=cos(phi)
      !sincos(4)=sin(phi)

!lambda = +1
      if(lambda.eq.1) then
         POL(1)= 0d0
         POL(2)= (sincos(3)*sincos(1)-(0d0,1d0)*sincos(4))/dsqrt(2d0)
         POL(3)= (sincos(4)*sincos(1)+(0d0,1d0)*sincos(3))/dsqrt(2d0)
         POL(4)= -sincos(2)/dsqrt(2d0)
!lambda = -1
      else if(lambda.eq.-1) then
         POL(1)= 0d0
         POL(2)= (sincos(3)*sincos(1)+(0d0,1d0)*sincos(4))/dsqrt(2d0)
         POL(3)= (sincos(4)*sincos(1)-(0d0,1d0)*sincos(3))/dsqrt(2d0)
         POL(4)= -sincos(2)/dsqrt(2d0)
!lambda = 0 (z)
      else if(lambda.eq.0) then
!|3-momentum|
         abs3p = dsqrt(p(2)**2+p(3)**2+p(4)**2)
!invariant mass
         inv_mass= dsqrt(p(1)**2-abs3p**2)
         POL(1)= abs3p/inv_mass
         POL(2)= sincos(3)*sincos(2)*p(1)/inv_mass
         POL(3)= sincos(4)*sincos(2)*p(1)/inv_mass
         POL(4)= sincos(1)*p(1)/inv_mass
!lambda = 2 (vec{p-hat})
      else if(lambda.eq.2) then
         POL(1)= 0d0
         POL(2)= sincos(3)*sincos(2)
         POL(3)= sincos(4)*sincos(2)
         POL(4)= sincos(1)
!lambda = -2 (-vec{p-hat})
      else if(lambda.eq.-2) then
         POL(1)= 0d0
         POL(2)= -sincos(3)*sincos(2)
         POL(3)= -sincos(4)*sincos(2)
         POL(4)= -sincos(1)
      endif

      return
      END subroutine POLARIZATION_SINGLE

!POLARIZATIONA.F
!VERSION 20130529

!in POL(lambda,^mu) returns the photon polarization vectors for given
!4-momentum p
      subroutine POLARIZATIONA(p, POL)

      implicit none
      real(8) p(4)
      complex(8) POL(2,4)

      call POLARIZATION_SINGLE(p,+1,POL(1,:))
      call POLARIZATION_SINGLE(p,-1,POL(2,:))

      return
      END subroutine POLARIZATIONA

!POLARIZATIONX.F
!VERSION 20130529

!in POL(lambda,^mu) returns the polarization vectors for given
!4-momentum p, with POL(L, ^mu) = (0, 0, 0, 1) when p is along the
!Z direction.
      subroutine POLARIZATIONX(p, POL)

      implicit none
      real(8) p(4)
      complex(8) POL(3,4)

      call POLARIZATION_SINGLE(p,+1,POL(1,:))
      call POLARIZATION_SINGLE(p,-1,POL(2,:))
      call POLARIZATION_SINGLE(p,+2,POL(3,:))

      return
      END subroutine POLARIZATIONX
!PROPAGATOR.F
!VERSION 20130522
!
!PROPAGATOR() returns the generi!complex-valued propagator
!without tensor structure (numerator), given mass, invariant mass
!and width.
      complex(8) function PROPAGATOR(inv_mass, mass, width)
      implicit none

      real(8) inv_mass, mass, width

!not assuming auto-conversion
!     PROPAGATOR = (0d0,1d0)/(dcmplx(inv_mass**2,0d0)
!    &            -dcmplx(mass**2,0d0)+
!    &             (0d0,1d0)*dcmplx(mass,0d0)*dcmplx(width,0d0))

!assuming auto-conversion. works with gfortran
      !print *,"Called propagator with",inv_mass,mass,width
      if (mass.ge.0d0) then
         PROPAGATOR = ci / ( inv_mass**2 - mass**2 + ci*mass*width )
      else
         PROPAGATOR = ci / ( inv_mass**2 )
      endif
!     print *, PROPAGATOR

      return
      END function PROPAGATOR

!VVP.F
!VERSION 20130524

!in epp(_mu, _nu) returns the
      subroutine VVP(p1,p2,epp)

      implicit none
      real(8) p1(4), p2(4)
      complex(8) epp(4,4)
!      real(8) ANTISYMMETRIC
      integer i,j,k,l

!      external ANTISYMMETRIC

      do i=1,4
      do j=1,4
        epp(i,j)=0d0
      enddo
      enddo

      do i=1,4
      do j=1,4
      do k=1,4
      do l=1,4
          epp(i,j)=epp(i,j)+ANTISYMMETRIC(i,j,k,l)*p1(k)*p2(l)
      enddo
      enddo
      enddo
      enddo

      return
      END subroutine VVP

!VVS1.F
!VERSION 20130524

!in g_mu_nu(_mu, _nu) returns the VVS1 4*4 tensor.
      subroutine VVS1(g_mu_nu)

      implicit none
      complex(8) g_mu_nu(4,4)

      g_mu_nu = 0d0

      g_mu_nu(1,1) =  1d0
      g_mu_nu(2,2) = -1d0
      g_mu_nu(3,3) = -1d0
      g_mu_nu(4,4) = -1d0

      return
      END subroutine VVS1

!VVS2.F
!VERSION 20130524

!in pp(_mu, _nu) returns the tensor of p1_mu * p2_nu.
      subroutine VVS2(p1,p2,pp)

      implicit none
      real(8) p1(4), p2(4)
      complex(8) pp(4,4)
      integer mu, nu

      do mu=1,4
        do nu=1,4
          pp(mu,nu)=p1(mu)*p2(nu)
          if( ( (mu.ne.1).and.(nu.eq.1) ).or. &
              ( (mu.eq.1).and.(nu.ne.1) ) )then
            pp(mu,nu)=-pp(mu,nu)
          endif
        enddo
      enddo

      return
      END subroutine VVS2

  !-- generic functions below
  !- MCFM spinors in non-MCFM momentum convention
      subroutine spinoru2(n,p,za,zb,s)
       implicit none
       integer, intent(in) :: n
       real(8), intent(in) :: p(4,n)
       complex(8), intent(out) :: za(n,n), zb(n,n)
       real(8), intent(out) :: s(n,n)
       integer :: i,j
       complex(8) :: c23(n), f(n)
       real(8) :: rt(n)

       !---if one of the vectors happens to be zero this routine fails.
       do j=1,N
          za(j,j)=czero
          zb(j,j)=za(j,j)

          !-----positive energy case
          if (p(1,j) .gt. zero) then
             rt(j)=sqrt(abs(p(2,j)+p(1,j)))
             c23(j)=dcmplx(p(4,j),-p(3,j))
             f(j)=(one,zero)
          else
          !-----negative energy case
             rt(j)=sqrt(abs(-p(1,j)-p(2,j)))
             c23(j)=dcmplx(-p(4,j),p(3,j))
             f(j)=ci
          endif
       enddo

       do i=2,N

        do j=1,i-1
             s(i,j)=two*scr(p(:,i),p(:,j))
             za(i,j)=f(i)*f(j)  * ( c23(i)*dcmplx(rt(j)/(rt(i)+1d-16))-c23(j)*dcmplx(rt(i)/(rt(j)+1d-16)) )

             if (abs(s(i,j)).lt.1d-5) then
                zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
             else
                zb(i,j)=-dcmplx(s(i,j))/(za(i,j)+1d-16)
             endif

             za(j,i)=-za(i,j)
             zb(j,i)=-zb(i,j)
             s(j,i)=s(i,j)

          enddo

       enddo

       return

      end subroutine spinoru2


!!!!!!!!!!!!!!
! XCHECK MEs !
!!!!!!!!!!!!!!

!MATRIXELEMENT02.F
!VERSION 20160924
!complex(8) function MATRIXELEMENT02(p,mass,helicity,id)
!use ModHiggs
!implicit none
!real(8), intent(in) :: p(1:4,1:9)
!real(8), intent(in) :: mass(3:5,1:2)
!real(8), intent(in) :: helicity(9)
!integer, intent(in) :: id(9)

!complex(dp) :: A_VV(1:4), propH
!integer :: MY_IDUP(6:9),i3,i4,j,VVMode
!real(dp) :: pUsed(4,6)
!integer :: ordering(1:4),ordering_swap(1:4)

!! Set the ids and momenta
!MY_IDUP(:)=Not_a_particle_
!pUsed(:,:)=zero
!do j=1,4
!   if(id(j) .ne. Not_a_particle_) then
!      if(j.le.2) then
!         MY_IDUP(j+5) = -convertLHEreverse(id(j))
!         pUsed(:,j+2) = -p(:,j)
!      else
!         MY_IDUP(j+5) = convertLHEreverse(id(j+3))
!         pUsed(:,j+2) = p(:,j+3)
!      endif
!   endif
!   pUsed(:,1) = pUsed(:,1) + pUsed(:,j+2)
!enddo
!call getDecay_VVMode_Ordering(MY_IDUP, VVMode,ordering,ordering_swap)

!! Set the helicities
!if(id(1)*helicity(1).gt.0) then
!   i3=2
!else
!   i3=1
!endif
!if(id(6)*helicity(6).gt.0) then
!   i4=2
!else
!   i4=1
!endif
!if(ordering(1).eq.5 .or. ordering(2).eq.5) then
!   call swap(i3,i4)
!endif

!A_VV(:) = 0d0
!PROPH = PROPAGATOR(dsqrt(scr(p(:,5),p(:,5))),mass(5,1),mass(5,2))
!call calcHelAmp2(ordering,VVMode,MY_IDUP,pUsed,i3,i4,A_VV(1))
!if( (VVMode.eq.ZZMode) .and. includeGammaStar ) then
!    call calcHelAmp2(ordering,ZgsMode,MY_IDUP,pUsed,i3,i4,A_VV(2))
!    call calcHelAmp2(ordering,gsZMode,MY_IDUP,pUsed,i3,i4,A_VV(3))
!    call calcHelAmp2(ordering,gsgsMode,MY_IDUP,pUsed,i3,i4,A_VV(4))
!elseif( VVMode.eq.ZgMode .and. includeGammaStar ) then
!    call calcHelAmp2(ordering,gsgMode,MY_IDUP,pUsed,i3,i4,A_VV(2))
!endif
!A_VV(:) = -A_VV(:) * propH

!MATRIXELEMENT02 = A_VV(1)+A_VV(2)+A_VV(3)+A_VV(4)
!return
!END function

!subroutine getDecay_VVMode_Ordering(MY_IDUP, VVMode,ordering,ordering_swap)
!   implicit none
!   integer, intent(in) :: MY_IDUP(6:9)
!   integer, intent(out) :: VVMode,ordering(1:4),ordering_swap(1:4)
!   integer :: idV(1:2)

!   ordering=(/3,4,5,6/)
!   idV(1)=CoupledVertex(MY_IDUP(6:7),-1)
!   idV(2)=CoupledVertex(MY_IDUP(8:9),-1)
!   if(MY_IDUP(6).eq.Pho_ .or. MY_IDUP(7).eq.Pho_) idV(1)=Pho_
!   if(MY_IDUP(8).eq.Pho_ .or. MY_IDUP(9).eq.Pho_) idV(2)=Pho_
!   if(convertLHE(MY_IDUP(6)).lt.0 .or. MY_IDUP(6).eq.Not_a_particle_) then
!      call swap(ordering(1),ordering(2))
!   endif
!   if(convertLHE(MY_IDUP(8)).lt.0 .or. MY_IDUP(8).eq.Not_a_particle_) then
!      call swap(ordering(3),ordering(4))
!   endif
!   if( &
!         (idV(1).eq.Wm_ .and. idV(2).eq.Wp_) .or. &
!         (idV(2).eq.Z0_ .and. idV(1).eq.Pho_) &
!     ) then
!      call swap(ordering(1),ordering(3))
!      call swap(ordering(2),ordering(4))
!      call swap(idV(1),idV(2))
!   endif
!   ordering_swap(:)=ordering(:)
!   call swap(ordering_swap(1),ordering_swap(3))

!   if(idV(1).eq.Z0_ .and. idV(2).eq.Z0_) then
!      VVMode=ZZMode
!   elseif(idV(1).eq.Z0_ .and. idV(2).eq.Pho_) then
!      VVMode=ZgMode
!   elseif(idV(1).eq.Pho_ .and. idV(2).eq.Pho_) then
!      VVMode=ggMode
!   elseif(idV(1).eq.Wp_ .and. idV(2).eq.Wm_) then
!      VVMode=WWMode
!   else
!      print *,"idV=",idV
!      call Error("Unsupported decay Modes")
!   endif
!   return
!end subroutine


SUBROUTINE MATRIXELEMENT1(p,FermFlav,UnPolSqAmp)
use ModParameters
use ModMisc
implicit none
complex(8) :: SME(1:3,-1:+1,-1:+1),HelAmp
real(8) :: p(1:4,1:9),UnpolSqAmp,PreFac,IZis(-1:+1)
real(8) :: qsq_V1,qsq_V2,qsq_V1V2,qsq_H
complex(8) ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
complex(8) :: a1HVV,a2HVV,a3HVV,Prop
integer :: ishel,fshel,FermFlav(1:6)! 12:IS, 34:ZDK, 56:HDK

    ! q1 qbar2 --> 3 --> 45 --> Z4-->f6 fbar7 + H5-->89
    call getSME(p,FermFlav,SME)
    if( H_DK ) call Error("Higgs decay not implemented")

    Prop = (0d0,1d0)/(((p(1:4,4)+p(1:4,5)).dot.(p(1:4,4)+p(1:4,5))) -M_V**2 + (0d0,1d0)*M_V*Ga_V )
    PreFac = 4d0*Pi*alpha_QED/4d0/sitW**2/(1d0-sitW**2)      ! gets squared below

    ! initial state couplings
    if( IsAWDecay(DecayMode1) ) then
        IZis(+1) = bR   *CKM(FermFlav(1),FermFlav(2))
        IZis(-1) = bL   *CKM(FermFlav(1),FermFlav(2))
    elseif( IsAZDecay(DecayMode1) .and. (abs(FermFlav(1)).eq.2 .or. abs(FermFlav(1)).eq.4) ) then
        IZis(+1) = aR_QUp
        IZis(-1) = aL_QUp
    elseif( IsAZDecay(DecayMode1) .and. (abs(FermFlav(1)).eq.1 .or. abs(FermFlav(1)).eq.3 .or. abs(FermFlav(1)).eq.5) ) then
        IZis(+1) = aR_QDn
        IZis(-1) = aL_QDn
    endif


    ! anomalous HVV couplings
    qsq_V1  =  p(1:4,3).dot.p(1:4,3)
    qsq_V2  =  p(1:4,4).dot.p(1:4,4)
    qsq_V1V2=-(p(1:4,3).dot.p(1:4,4))
    qsq_H   =  p(1:4,5).dot.p(1:4,5)

    ghz1_dyn = HVVSpinZeroDynamicCoupling(1,qsq_V1,qsq_V2,qsq_H)
    ghz2_dyn = HVVSpinZeroDynamicCoupling(2,qsq_V1,qsq_V2,qsq_H)
    ghz3_dyn = HVVSpinZeroDynamicCoupling(3,qsq_V1,qsq_V2,qsq_H)
    ghz4_dyn = HVVSpinZeroDynamicCoupling(4,qsq_V1,qsq_V2,qsq_H)

    a1HVV = ghz1_dyn*M_V**2 + qsq_V1V2*( 2d0*ghz2_dyn + ghz3_dyn*qsq_V1V2/Lambda**2 )
    a2HVV =-2d0*ghz2_dyn - ghz3_dyn*qsq_V1V2/Lambda**2
    a3HVV =-2d0*ghz4_dyn

    UnPolSqAmp = 0d0
    do ishel=-1,+1,2
    do fshel=-1,+1,2
          HelAmp =  a3HVV * (  - IZis(ishel)*SME(3,ishel,fshel)  )          &
                  + a2HVV * (  - IZis(ishel)*SME(2,ishel,fshel)  )          &
                  + a1HVV * (  + IZis(ishel)*SME(1,ishel,fshel)  )
          HelAmp = HelAmp * PreFac/vev * Prop
          UnPolSqAmp = UnPolSqAmp + dreal( HelAmp*dconjg(HelAmp) )
    enddo
    enddo
    UnPolSqAmp = UnPolSqAmp * CF

RETURN
END SUBROUTINE

SUBROUTINE getSME(p,FermFlav,SME)
use ModParameters
use ModMisc
implicit none
complex(8) :: SME(1:3,-1:+1,-1:+1)
real(8) :: sprod(9,9),p(1:4,1:9),IZfs(-1:+1)
complex(8) :: za(9,9), zb(9,9),Prop
integer :: FermFlav(1:6)

    call spinoru2(9,(/-p(1:4,1),-p(1:4,2),-p(1:4,1)-p(1:4,2),p(1:4,6)+p(1:4,7),p(1:4,8)+p(1:4,9),p(1:4,6),p(1:4,7),p(1:4,8),p(1:4,9)/),za,zb,sprod)


    ! Z-final state couplings
    if( IsAWDecay(DecayMode1) ) then
        IZfs(+1) = bR   *CKM(FermFlav(3),FermFlav(4))
        IZfs(-1) = bL   *CKM(FermFlav(3),FermFlav(4))
    elseif( abs(FermFlav(3)).eq.11 .or. abs(FermFlav(3)).eq.13 .or. abs(FermFlav(3)).eq.15) then
         IZfs(-1)=aL_lep    * dsqrt(scale_alpha_Z_ll)
         IZfs(+1)=aR_lep    * dsqrt(scale_alpha_Z_ll)
    elseif( abs(FermFlav(3)).eq.12 .or. abs(FermFlav(3)).eq.14 .or. abs(FermFlav(3)).eq.16 ) then
         IZfs(-1)=aL_neu    * dsqrt(scale_alpha_Z_nn)
         IZfs(+1)=aR_neu    * dsqrt(scale_alpha_Z_nn)
    elseif( abs(FermFlav(3)).eq.2 .or. abs(FermFlav(3)).eq.4 ) then
         IZfs(-1)=aL_QUp    * dsqrt(scale_alpha_Z_uu)
         IZfs(+1)=aR_QUp    * dsqrt(scale_alpha_Z_uu)
    elseif( abs(FermFlav(3)).eq.1 .or. abs(FermFlav(3)).eq.3 .or. abs(FermFlav(3)).eq.5 ) then
         IZfs(-1)=aL_QDn    * dsqrt(scale_alpha_Z_dd)
         IZfs(+1)=aR_QDn    * dsqrt(scale_alpha_Z_dd)
    else
         call Error("Wrong flavor in getSME",FermFlav(3))
    endif

    SME(1,+1,+1) =  -2*IZfs(1)*za(1,7)*zb(2,6)
    SME(2,+1,+1) =  IZfs(1)*(za(1,6)*zb(2,6) + za(1,7)*zb(2,7))*(za(7,8)*zb(6,8) + za(7,9)*zb(6,9))
    SME(3,+1,+1) = cI*IZfs(1)*(za(1,7)*(za(7,8)*zb(2,8) + za(7,9)*zb(2,9))*zb(6,7) + za(6,7)*zb(2,6)*(za(1,8)*zb(6,8) + za(1,9)*zb(6,9)))
    SME(1,+1,-1) =  -2*IZfs(-1)*za(1,6)*zb(2,7)
    SME(2,+1,-1) =  IZfs(-1)*(za(1,6)*zb(2,6) + za(1,7)*zb(2,7))*(za(6,8)*zb(7,8) + za(6,9)*zb(7,9))
    SME(3,+1,-1) = cI*IZfs(-1)*(-(za(1,7)*zb(2,7)*(za(6,8)*zb(7,8) + za(6,9)*zb(7,9))) + za(1,6)*(za(6,8)*(-(zb(2,7)*zb(6,8)) + zb(2,6)*zb(7,8)) + zb(2,7)*(za(7,8)*zb(7,8) + za(7,9)*zb(7,9)) + za(6,9)*(-(zb(2,7)*zb(6,9)) + zb(2,6)*zb(7,9))))
    SME(1,-1,+1) = -2*IZfs(1)*za(2,7)*zb(1,6)
    SME(2,-1,+1) = IZfs(1)*(za(2,6)*zb(1,6) + za(2,7)*zb(1,7))*(za(7,8)*zb(6,8) + za(7,9)*zb(6,9))
    SME(3,-1,+1) = cI*IZfs(1)*(za(2,7)*(za(7,8)*zb(1,8) + za(7,9)*zb(1,9))*zb(6,7) + za(6,7)*zb(1,6)*(za(2,8)*zb(6,8) + za(2,9)*zb(6,9)))
    SME(1,-1,-1) = -2*IZfs(-1)*za(2,6)*zb(1,7)
    SME(2,-1,-1) =  IZfs(-1)*(za(2,6)*zb(1,6) + za(2,7)*zb(1,7))*(za(6,8)*zb(7,8) + za(6,9)*zb(7,9))
    SME(3,-1,-1) =   -(cI*IZfs(-1)*(za(2,6)*(za(6,8)*zb(1,8) + za(6,9)*zb(1,9))*zb(6,7) + za(6,7)*zb(1,7)*(za(2,8)*zb(7,8) + za(2,9)*zb(7,9))))

    Prop = (0d0,1d0)/(2*(p(1:4,6).dot.p(1:4,7)) - M_V**2 + (0d0,1d0)*M_V*Ga_V )
    SME(:,:,:) = SME(:,:,:) * Prop

RETURN
END SUBROUTINE







      ! p: g1, g2, Zl1, Zl2, H
      SUBROUTINE EvalUnpolAmpSq_gg_VH(p,UnPolSqAmp)
      use ModParameters
      use ModMisc
#if useCollier==1
      use COLLIER
#endif
      implicit none
      real(8) :: p(1:4,1:4),UnPolSqAmp
      complex(8) :: SP(1:7),SI(0:16),LC(1:15)
      complex(8) :: qlI1,qlI2,qlI3,qlI4,testYA,testMA
      complex(8) :: ep1(1:4),ep2(1:4),cep3(1:4),p1(1:4),p2(1:4),p3(1:4)
      complex(8) :: MasslessTri,MassiveTri,MassiveBox
      real(8) :: shat,that,uhat,MT2,MH2,MZ2,MV2
      real(8) :: MuRen2,PreFac,PreFac2,IZ(3:4,9:9),smT,y,mTrans
      integer :: eps,h1,h2,h3,i,j
      complex(8) :: VVHg1,VVHg2,VVHg3
      complex(8) :: cmzsq, rat,cix,newsum_MCFM
      complex(8) :: calc_MassiveBox,calc_MassiveBox_QP,calc_MassiveTensorBox,newsum
      real(8),parameter :: E=dexp(1d0), ICol_sq=8d0
      integer,parameter :: rank=3
      double complex :: DD1(0:rank/2,0:rank,0:rank,0:rank),DD1uv(0:rank/2,0:rank,0:rank,0:rank)
      double complex :: DD2(0:rank/2,0:rank,0:rank,0:rank),DD2uv(0:rank/2,0:rank,0:rank,0:rank)
      double complex :: DD3(0:rank/2,0:rank,0:rank,0:rank),DD3uv(0:rank/2,0:rank,0:rank,0:rank)
      double complex :: DD0(1:3)
      double precision :: Derr(0:rank)
      real(8) :: sprod(1:4,1:4),p5hat(1:4),p6hat(1:4)
      complex(8) :: za(1:4,1:4),zb(1:4,1:4),a,b,c,d,ee,f,v,w,x,yy,z,o,ampMCFM(-1:1,-1:1,-1:1)

!       prefactors and kinematics
        PreFac  = Pi**2 *      (4d0*pi*alpha_QED) * (4d0*pi*alphas) /sitW/M_W     *1d0/(2d0*pi)**4
        PreFac2 = Pi**2 * dsqrt(4d0*pi*alpha_QED) * (4d0*pi*alphas) *m_Top**2/vev *1d0/(2d0*pi)**4
        IZ(3,9) = 0.5d0*(aL_QUp-aR_QUp)/2d0/sitW/dsqrt(1d0-sitW**2) ! top
        IZ(4,9) = -IZ(3,9)                                          ! bottom

        p1(1:4) = p(1:4,1)
        p2(1:4) = p(1:4,2)
        p3(1:4) = p(1:4,3)+p(1:4,4)
        MZ2 = M_Z**2
        MH2 = M_Reso**2
        MT2 = M_Top**2
        MV2 = get_MInv2(dble(p3(1:4)))
        shat =+2d0*(p1.dot.p2)
        that =-2d0*(p1.dot.p3)+MV2
        uhat =-2d0*(p2.dot.p3)+MV2
        mTrans = dsqrt( get_pT(dble(p3(1:4)))**2 + MV2 )
        smT  = dsqrt(shat)*mTrans
        y    = Get_ETA(dble(p3(1:4)))

        call spinoru2(4,(/p(1:4,1),p(1:4,2),p(1:4,3),p(1:4,4)/),za,zb,sprod)
        a = za(1,2);  b = za(1,3);  c = za(1,4);  d = za(2,3);  ee = za(2,4);  f = za(3,4)
        v = zb(1,2);  w = zb(1,3);  x = zb(1,4);  yy = zb(2,3);  z = zb(2,4);  o = zb(3,4)


!       evaluate scalar integrals
        SI(0) = 1d0
        eps = 0
!         ! check 1/eps^2
!         print *, "checking 1/eps^2"
!         eps=-2
!         SI(0) = 0d0
        ! check 1/eps^1
!         print *, "checking 1/eps^1"
!         eps=-1
!         SI(0) = 0d0

        MuRen2 = Mu_Ren**2
!         print *, "checks"
!         print *, PreFac,IZ,MZ2,MH2,MT2
!         print *, shat,that,uhat
!         print *, p(1:4,1)
!         print *, p(1:4,2)
!         print *, p(1:4,3)
!         print *, p(1:4,4)
!         pause
#if useCollier==1
         call SetMuUV2_cll(MuRen2)
         call SetMuIR2_cll(MuRen2)
#endif


!         SI(1:16) = 0d0
#if useQCDLoopLib==1 || linkMELA==1
        SI(1) = qlI1(MT2,MuRen2,eps)
        SI(2) = qlI2(zero, MT2,MT2,MuRen2,eps)
        SI(3) = qlI2(MH2,  MT2,MT2,MuRen2,eps)
        SI(4) = qlI2(MV2,  MT2,MT2,MuRen2,eps)
        SI(5) = qlI2(shat, MT2,MT2,MuRen2,eps)
        SI(6) = qlI2(that, MT2,MT2,MuRen2,eps)
        SI(7) = qlI2(uhat, MT2,MT2,MuRen2,eps)
        SI(8) = qlI3(zero,zero,shat, MT2,MT2,MT2,MuRen2,eps)
        SI(9) = qlI3(zero,MV2,that, MT2,MT2,MT2,MuRen2,eps)
        SI(10)= qlI3(zero,MV2,uhat, MT2,MT2,MT2,MuRen2,eps)
        SI(11)= qlI3(zero,that,MH2, MT2,MT2,MT2,MuRen2,eps)
        SI(12)= qlI3(zero,uhat,MH2, MT2,MT2,MT2,MuRen2,eps)
        SI(13)= qlI3(MV2,shat,MH2,  MT2,MT2,MT2,MuRen2,eps)
        SI(14)= qlI4(zero,zero,MV2,MH2,shat,that, MT2,MT2,MT2,MT2,MuRen2,eps)
        SI(15)= qlI4(zero,zero,MV2,MH2,shat,uhat, MT2,MT2,MT2,MT2,MuRen2,eps)
        SI(16)= qlI4(zero,MV2,zero,MH2,that,uhat, MT2,MT2,MT2,MT2,MuRen2,eps)
#else
    call Error("This process requires linkQCDLoop or linkMELA enabled in the makefile")
#endif


#if useCollier==1
      call D_cll(DD1,DD1uv,dcmplx((/0d0,MH2,MV2,0d0,that,shat/)),dcmplx((/MT2,MT2,MT2,MT2/)),rank,Derr)  
      call D_cll(DD2,DD2uv,dcmplx((/0d0,MV2,MH2,0d0,uhat,shat/)),dcmplx((/MT2,MT2,MT2,MT2/)),rank,Derr)  
      call D_cll(DD3,DD3uv,dcmplx((/MV2,0d0,MH2,0d0,uhat,that/)),dcmplx((/MT2,MT2,MT2,MT2/)),rank,Derr)      

      call D0_cll(DD0(1),dcmplx((/0d0,0d0,MV2,MH2,shat,that/)),dcmplx((/MT2,MT2,MT2,MT2/)))  
      call D0_cll(DD0(2),dcmplx((/0d0,0d0,MV2,MH2,shat,uhat/)),dcmplx((/MT2,MT2,MT2,MT2/)))  
      call D0_cll(DD0(3),dcmplx((/0d0,MV2,0d0,MH2,that,uhat/)),dcmplx((/MT2,MT2,MT2,MT2/)))
#else
      call Error("This process requires linkCollier enabled in the makefile")
#endif


             
 VVHg1 = 0.5d0*ghz1 + ghz2 * ( p3(1:4).dot.(p1(1:4)+p2(1:4)) )/MZ2
 VVHg2 = -ghz2/MZ2
 

! !      triangles in spinor helicity (for one helicity configuration)
!        newsum = (-4*a*(0d0,1d0)*MT2*VVHg1*(c*w + ee*yy)*SI(8))/v       
!        newsum = newsum  &
!               + (-2*a*(0d0,1d0)*MT2*VVHg2/MH2*(b**2*w*x + b*(d*o*v - 2*a*v*x + c*x**2 + 2*d*x*yy) + z*(a*(-2*d*v + f*x) + d*(2*c*x + d*yy + ee*z))))*SI(8)/v
!        newsum = newsum /( MV2-MZ2 + (0d0,1d0)*M_Z*Ga_Z )   ! * exp((0d0,1d0)*Pi*(+1d0))

       
       
!        print *, "TESTER",(MZ2 + a*v)/(-MH2 + MV2 + MZ2 + b*w + c*x + d*yy + ee*z)



! ! ! ! ! ! ! MCFM triangle  (remember to set alpha_s=constant)
! !------ right lepton amplitudes
!       amp(2,2,2)=ggHZ_pp_tri(1,2,4,3,za,zb,mt2)

!       call spinoru2(4,(/-p(1:4,1),-p(1:4,2),p(1:4,4),p(1:4,3)/),za,zb,sprod)
!       cmzsq=M_Z**2-(0d0,1d0)*M_Z*Ga_Z
! ! top tri
!       cix=  (2*(-((mt2*za(1,3)*zb(2,1)*zb(4,1))/za(1,2)) -  (mt2*za(2,3)*zb(2,1)*zb(4,2))/za(1,2)))/ (za(1,2)*za(3,4)*zb(2,1)*zb(4,3))
!       cix=cix+ (2*mt2*zb(2,1)**2*(za(1,3)*zb(4,1) + za(2,3)*zb(4,2)))/  (cmzsq*sprod(1,2)*sprod(3,4))
!       rat= (-2*(-((za(1,3)*zb(2,1)*zb(4,1))/za(1,2)) -   (za(2,3)*zb(2,1)*zb(4,2))/za(1,2)))/(za(1,2)*za(3,4)*zb(2,1)*zb(4,3))
!       rat=rat+ (-2*zb(2,1)**2* (za(1,3)*zb(4,1) + za(2,3)*zb(4,2)))/ (cmzsq*sprod(1,2)*sprod(3,4))
!       newsum_MCFM=( SI(8)*cix-rat/2d0 ) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z) * sprod(1,2)/dcmplx(sprod(1,2)-M_Z**2,0d0) 
! 
! ! bot tri
! mt2=0d0
! SI(8) = qlI3(zero,zero,shat, MT2,MT2,MT2,MuRen2,eps)
!       cix=  (2*(-((mt2*za(1,3)*zb(2,1)*zb(4,1))/za(1,2)) -  (mt2*za(2,3)*zb(2,1)*zb(4,2))/za(1,2)))/ (za(1,2)*za(3,4)*zb(2,1)*zb(4,3))
!       cix=cix+ (2*mt2*zb(2,1)**2*(za(1,3)*zb(4,1) + za(2,3)*zb(4,2)))/  (cmzsq*sprod(1,2)*sprod(3,4))
!       rat= (-2*(-((za(1,3)*zb(2,1)*zb(4,1))/za(1,2)) -   (za(2,3)*zb(2,1)*zb(4,2))/za(1,2)))/(za(1,2)*za(3,4)*zb(2,1)*zb(4,3))
!       rat=rat+ (-2*zb(2,1)**2* (za(1,3)*zb(4,1) + za(2,3)*zb(4,2)))/ (cmzsq*sprod(1,2)*sprod(3,4))
!       newsum_MCFM=newsum_MCFM - ( SI(8)*cix-rat/2d0 ) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z) * sprod(1,2)/dcmplx(sprod(1,2)-M_Z**2,0d0) 

! restoring       
! mt2 = m_top**2     
! SI(8) = qlI3(zero,zero,shat, MT2,MT2,MT2,MuRen2,eps)
! 
! 
! 
! 
! ! the boxes 
!       call spinoru2(4,(/-p(1:4,1),-p(1:4,2),p(1:4,4),p(1:4,3)/),za,zb,sprod)
! 
! !------ left handed lepton coupling (- sign from line reversal)
!       ampMCFM(-1,-1,+1)= ggHZ_pp_box(1,2,3,4,za,zb,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aL_lep
!       ampMCFM(+1,-1,+1)= ggHZ_mp_box(1,2,3,4,za,zb,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aL_lep
!       ampMCFM(-1,+1,+1)=-ggHZ_mp_box(1,2,4,3,zb,za,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aL_lep
!       ampMCFM(+1,+1,+1)=-ggHZ_pp_box(1,2,4,3,zb,za,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aL_lep
!       
! !------- right handed lepton coupling 
!       ampMCFM(-1,-1,-1)= ggHZ_pp_box(1,2,4,3,za,zb,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aR_lep
!       ampMCFM(+1,-1,-1)= ggHZ_mp_box(1,2,4,3,za,zb,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aR_lep
!       ampMCFM(-1,+1,-1)=-ggHZ_mp_box(1,2,3,4,zb,za,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aR_lep
!       ampMCFM(+1,+1,-1)=-ggHZ_pp_box(1,2,3,4,zb,za,sprod,mt2) *sprod(3,4)/dcmplx(sprod(3,4)-M_Z**2,M_Z*Ga_Z)  * aR_lep
!       
! ! ! ! ! ! END MCMF
            


!       helicity sum
        UnPolSqAmp = 0d0
        do h1=-1,+1, 2
        do h2=-1,+1, 2
        do h3=-1,+1, 2

!           polarization vectors
            ep1(1:4) = pol_gluon_incoming(p(1:4,1),h1)
            ep2(1:4) = pol_gluon_incoming(p(1:4,2),h2)

!             cep3(1:4) = pol_mass(p3(1:4),h3,outgoing=.true.)

            
            cep3(1:4)= pol_Zff_outgoing(p(1:4,3),p(1:4,4),h3)   * (-1d0)/( MV2-MZ2 + (0d0,1d0)*M_Z*Ga_Z )
            if( h3.lt.0 ) then
               cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aR_lep
            else
               cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aL_lep
            endif



!           scalar products
            SP(1) = (cep3.dot.ep1)
            SP(2) = (cep3.dot.ep2)
            SP(3) = (cep3.dot.p1)
            SP(4) = (cep3.dot.p2)
            SP(5) = ( ep1.dot.ep2)
            SP(6) = ( ep1.dot.p3)
            SP(7) = ( ep2.dot.p3)

!           Levi-Civita tensors
            LC(1) = LeviCiv(p1,ep1,ep2,cep3)
            LC(2) = LeviCiv(p1,p2,ep1,cep3)
            LC(3) = LeviCiv(p1,p2,ep1,ep2)
            LC(4) = LeviCiv(p1,p2,ep2,cep3)
            LC(5) = LeviCiv(p1,p2,p3,cep3)
            LC(6) = LeviCiv(p1,p2,p3,ep1)
            LC(7) = LeviCiv(p1,p2,p3,ep2)
            LC(8) = LeviCiv(p1,p3,ep1,cep3)
            LC(9) = LeviCiv(p1,p3,ep1,ep2)
            LC(10)= LeviCiv(p1,p3,ep2,cep3)
            LC(11)= LeviCiv(p2,ep1,ep2,cep3)
            LC(12)= LeviCiv(p2,p3,ep1,cep3)
            LC(13)= LeviCiv(p2,p3,ep1,ep2)
            LC(14)= LeviCiv(p2,p3,ep2,cep3)
            LC(15)= LeviCiv(p3,ep1,ep2,cep3)

!             shat*(LC(1) - LC(11)) - 2*LC(3)*(SP(3) + SP(4))  == 00 

print *, "CHECK: for some reasone the original VVHg2,VVHg3 couplings were not divided by MH2"
print *, "CHECK in comparison with Yaofu/Mathematica because VVHg2 is already divided above."

!           Loop amplitudes
            MasslessTri = (2*MZ2*SI(0)*(-4*VVHg1*LC(3)*(SP(3) + SP(4)) + 4*smT*VVHg2/MH2*Cosh(y)*LC(3)*(SP(3) + SP(4))  & 
                            + shat*(-(VVHg1*(LC(1) + 5*LC(11))) + VVHg2/MH2*(LC(9) + 5*LC(13))*(SP(3) + SP(4)))         &
                            + 3*shat**2*VVHg3/MH2*(SP(2)*SP(6) - SP(1)*SP(7))))/(3*shat*(-MH2 + MV2 + MZ2 - 2*smT*Cosh(y)))
!
            MassiveTri = (-2*(2*smT*VVHg2/MH2*LC(3)*(MZ2*(shat*SI(0) - SI(1) + MT2*(SI(0) - 2*SI(2) + 3*SI(5))) +   &
                           3*MT2*(MZ2 - shat)*shat*SI(8))*(SP(3) + SP(4)) + 2*E**(2*y)*smT*VVHg2/MH2*LC(3)*   &
                           (MZ2*(shat*SI(0) - SI(1) + MT2*(SI(0) - 2*SI(2) + 3*SI(5))) + 3*MT2*(MZ2 - shat)*shat*SI(8))*(SP(3) + SP(4)) +    &
                           E**y*(-2*MT2*(-6*shat**2*(VVHg1 + shat*VVHg2/MH2)*LC(3)*SI(8)*(SP(3) + SP(4)) + MZ2*(2*VVHg1*LC(3)*(SI(0) - 2*SI(2) + 3*SI(5))*   &
                           (SP(3) + SP(4)) + 6*shat**2*VVHg2/MH2*LC(3)*SI(8)*(SP(3) + SP(4)) + shat*(VVHg2/MH2*(LC(9) - LC(13))*(SI(0) - 2*SI(2) + 3*SI(5))*   &
                           (SP(3) + SP(4)) + VVHg1*(-((LC(1) - LC(11))*(SI(0) - 2*SI(2) + 3*SI(5))) +    &
                           6*LC(3)*SI(8)*(SP(3) + SP(4)))))) +  MZ2*(4*VVHg1*LC(3)*SI(1)*(SP(3) + SP(4)) +    &
                           shat**2*SI(0)*(-(VVHg1*(LC(1) + 5*LC(11))) +  VVHg2/MH2*(LC(9) + 5*LC(13))*(SP(3) + SP(4))) -    &
                           2*shat*(-(VVHg2/MH2*(LC(9) - LC(13))*SI(1)*(SP(3) + SP(4))) + VVHg1*((LC(1) - LC(11))*SI(1) +    &
                           2*LC(3)*SI(0)*(SP(3) + SP(4)))) + 3*shat**3*VVHg3/MH2*SI(0)*(SP(2)*SP(6) - SP(1)*SP(7))))))/(3*E**y*shat**2*(MH2 - MV2 - MZ2 + 2*smT*Cosh(y)))


!             MassiveBox = calc_MassiveBox(MV2,MH2,MT2,shat,smT,y,LC,SP,SI,kappa,kappa_tilde)
!             MassiveBox = calc_MassiveBox_QP(MV2,MH2,MT2,shat,smT,y,LC,SP,SI,kappa,kappa_tilde)
            MassiveBox = calc_MassiveTensorBox(MV2,MT2,MH2,shat,that,uhat,smT,y,LC,SP,DD0,DD1,DD2,DD3,kappa,kappa_tilde)

            MasslessTri = MasslessTri * IZ(4,9)* PreFac 
            MassiveTri  = MassiveTri  * IZ(3,9)* PreFac 
            MassiveBox  = MassiveBox  * IZ(3,9)* PreFac2

            UnPolSqAmp = UnPolSqAmp + ICol_sq * cdabs( MasslessTri*1 + MassiveTri*1 + MassiveBox*1 )**2



!             print *, h1,h2,h3
!             print *, "MasslessTri",MasslessTri
!             print *, "MassiveTri",MassiveTri
!             print *, "MassiveBox",MassiveBox
!             pause

        enddo
        enddo
        enddo


      RETURN
      END SUBROUTINE






FUNCTION TreeAmp_qqb_ZH(MomExt,Heli)! ordering: 1:in 2:in 3:H 4:l1 5:l2
use ModParameters
use ModMisc
implicit none
integer,parameter:: minus=0, plus=1, dummy=0
real(8) :: MomExt(1:4,1:6),MomZst1(1:4),MomZst2(1:4),PreFac,MomH(1:4)
integer :: Heli(1:4)
complex(8) ::TreeAmp_qqb_ZH(1:6),TreeAmp(minus:plus),VVHg1,VVHg2,VVHg3
complex(8) :: Spi(1:4),BarSpi(1:4,minus:plus),ZPol1(1:4,minus:plus),ZPol2(1:4),cep3(1:4),ZProp1,ZProp2
TreeAmp_qqb_ZH(1:6) = (0d0,0d0)


    Spi(1:4) = v_spinor(MomExt(1:4,1),-Heli(1))             !  incomming quark= u(p,lambda) = v(p,-lambda)
    BarSpi(1:4,dummy) = ubar_spinor(MomExt(1:4,2),-Heli(2)) !  incomming anti-quark= vbar(p,lambda) = ubar(p,-lambda)
    BarSpi(1:4,plus)  = Chir_Weyl(plus,BarSpi(1:4,dummy))
    BarSpi(1:4,minus) = Chir_Weyl(minus,BarSpi(1:4,dummy))
    ZPol1(1:4,plus)  = vbqq_Weyl(BarSpi(1:4,plus),Spi(1:4))  
    ZPol1(1:4,minus) = vbqq_Weyl(BarSpi(1:4,minus),Spi(1:4)) 



    MomZst1(1:4) = MomExt(1:4,1)+MomExt(1:4,2)
    MomZst2(1:4) = MomExt(1:4,5)+MomExt(1:4,6)
    MomH(1:4)    = MomExt(1:4,3)+MomExt(1:4,4)    
    ZProp1 = -(0d0,1d0)/((MomZst1(1:4).dot.MomZst1(1:4))-M_Z**2 +(0d0,1d0)*Ga_Z*M_Z)
    ZProp2 = -(0d0,1d0)/((MomZst2(1:4).dot.MomZst2(1:4))-M_Z**2 +(0d0,1d0)*Ga_Z*M_Z)
    

    cep3(1:4)= pol_Zff_outgoing(MomExt(1:4,5),MomExt(1:4,6),Heli(3)) * ZProp2
    if( Heli(3).lt.0 ) then
       cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aR_lep
    else
       cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aL_lep
    endif    

    PreFac = 2d0*m_v**2/vev 
    VVHg1 = PreFac*( 0.5d0*ghz1 + ghz2 * ( MomZst2.dot.MomZst1 )/M_Z**2 )
    VVHg2 = PreFac*( -ghz2/M_Z**2 )
    VVHG3 = PreFac*( -ghz3/M_Z**2 )
    
    
    TreeAmp(plus)  =  VVHg1 * (cep3.dot.ZPol1(:,plus))  + VVHg2*(MomH(1:4).dot.cep3) * (ZPol1(:,plus).dot.MomH(1:4))  + VVHg3 * LeviCiv(ZPol1(:,plus), cep3,dcmplx(MomZst2),dcmplx(MomH(1:4)))
    TreeAmp(minus) =  VVHg1 * (cep3.dot.ZPol1(:,minus)) + VVHg2*(MomH(1:4).dot.cep3) * (ZPol1(:,minus).dot.MomH(1:4)) + VVHg3 * LeviCiv(ZPol1(:,minus),cep3,dcmplx(MomZst2),dcmplx(MomH(1:4)))
    
    TreeAmp_qqb_ZH(up_)  = (TreeAmp(plus)*aR_QUp + TreeAmp(minus)*aL_QUp) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqb_ZH(dn_)  = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqb_ZH(chm_) = (TreeAmp(plus)*aR_QUp + TreeAmp(minus)*aL_QUp) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqb_ZH(str_) = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqb_ZH(bot_) = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1

RETURN
END FUNCTION






      

FUNCTION TreeAmp_qqbg_ZH(MomExt,Heli)! gluon emission from initial state ! ordering: 1:in 2:in 3:H 4:l1 5:l2 6: gluon
use ModParameters
use ModMisc
implicit none
integer,parameter:: minus=-1, plus=1, dummy=0
real(8) :: MomExt(1:4,1:7)
integer :: Heli(1:5)
complex(8) ::TreeAmp_qqbg_ZH(1:6),TreeAmp(minus:plus),VVHg1,VVHg2,VVHg3,cep3(1:4)
complex(8) :: Spi(1:4),SpiAux(1:4,minus:plus),PolGlu(1:4),BarSpi(1:4),ZPol1(1:4,minus:plus),ZPol2(1:4),ZProp1,ZProp2,PropMom(1:4)
real(8) :: MomZst1(1:4),MomZst2(1:4),MomH(1:4),PreFac
TreeAmp_qqbg_ZH(:) = (0d0,0d0)


    Spi(1:4) = v_spinor(MomExt(1:4,1),-Heli(1))             !  incomming quark= u(p,lambda) = v(p,-lambda)
    BarSpi(1:4) = ubar_spinor(MomExt(1:4,2),-Heli(2))       !  incomming anti-quark= vbar(p,lambda) = ubar(p,-lambda)
    PolGlu(1:4) = pol_mless(dcmplx(MomExt(1:4,7)),Heli(5),outgoing=.true.)
! PolGlu(1:4) = MomExt(1:4,7) ; print *, "gauge invariance check"
    
    MomZst1(1:4) = MomExt(1:4,1)+MomExt(1:4,2)-MomExt(1:4,7) 
    MomZst2(1:4) = MomExt(1:4,5)+MomExt(1:4,6)
    MomH(1:4)    = MomExt(1:4,3)+MomExt(1:4,4)

!   diagram 1: gluon emission from quark
    SpiAux(1:4,dummy) = spi2_Weyl(PolGlu(1:4),Spi(1:4))* cI * dsqrt(4d0*Pi*alphas)
    PropMom(1:4) = MomExt(1:4,1)-MomExt(1:4,7)
    SpiAux(1:4,dummy) = spi2_Weyl(PropMom(1:4),SpiAux(1:4,dummy)) * cI/(PropMom(1:4).dot.PropMom(1:4))
    SpiAux(1:4,plus)  = Chir_Weyl(plus, BarSpi(1:4))
    SpiAux(1:4,minus) = Chir_Weyl(minus,BarSpi(1:4))
    ZPol1(1:4,plus)   = vbqq_Weyl(SpiAux(1:4,plus), SpiAux(1:4,dummy))  
    ZPol1(1:4,minus)  = vbqq_Weyl(SpiAux(1:4,minus),SpiAux(1:4,dummy)) 

!   diagram 2: gluon emisson from anti-quark
    SpiAux(1:4,dummy) = spb2_Weyl(BarSpi(1:4),PolGlu(1:4))* cI * dsqrt(4d0*Pi*alphas)
    PropMom(1:4)=-(MomExt(1:4,2)-MomExt(1:4,7))
    SpiAux(1:4,dummy) = spb2_Weyl(SpiAux(1:4,dummy),PropMom(1:4)) * cI/(PropMom(1:4).dot.PropMom(1:4))
    SpiAux(1:4,plus)  = Chir_Weyl(plus, SpiAux(1:4,dummy))
    SpiAux(1:4,minus) = Chir_Weyl(minus,SpiAux(1:4,dummy))
    ZPol1(1:4,plus)   = ZPol1(1:4,plus)  + vbqq_Weyl(SpiAux(1:4,plus), Spi(1:4)) 
    ZPol1(1:4,minus)  = ZPol1(1:4,minus) + vbqq_Weyl(SpiAux(1:4,minus),Spi(1:4)) 
    
    ZProp1 = -(0d0,1d0)/((MomZst1(1:4).dot.MomZst1(1:4))-M_Z**2 +(0d0,1d0)*Ga_Z*M_Z)
    ZProp2 = -(0d0,1d0)/((MomZst2(1:4).dot.MomZst2(1:4))-M_Z**2 +(0d0,1d0)*Ga_Z*M_Z)
    
    cep3(1:4)= pol_Zff_outgoing(MomExt(1:4,5),MomExt(1:4,6),Heli(3)) * ZProp2
    if( Heli(3).lt.0 ) then
       cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aR_lep
    else
       cep3(1:4) = cep3(1:4) * dsqrt(couplZffsq) * aL_lep
    endif    
     
             
    PreFac = 2d0*m_v**2/vev 
    VVHg1 = PreFac*( 0.5d0*ghz1 + ghz2 * ( MomZst2.dot.MomZst1 )/M_Z**2 )
    VVHg2 = PreFac*( -ghz2/M_Z**2 )
    VVHG3 = PreFac*( -ghz3/M_Z**2 )
     
    
    TreeAmp(plus)  =  VVHg1 * (cep3.dot.ZPol1(:,plus))  + VVHg2*(MomH(1:4).dot.cep3) * (ZPol1(:,plus).dot.MomH(1:4))  + VVHg3 * LeviCiv(ZPol1(:,plus), cep3,dcmplx(MomZst2),dcmplx(MomH(1:4)))
    TreeAmp(minus) =  VVHg1 * (cep3.dot.ZPol1(:,minus)) + VVHg2*(MomH(1:4).dot.cep3) * (ZPol1(:,minus).dot.MomH(1:4)) + VVHg3 * LeviCiv(ZPol1(:,minus),cep3,dcmplx(MomZst2),dcmplx(MomH(1:4)))
    

    TreeAmp_qqbg_ZH(up_)  = (TreeAmp(plus)*aR_QUp + TreeAmp(minus)*aL_QUp) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqbg_ZH(dn_)  = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqbg_ZH(chm_) = (TreeAmp(plus)*aR_QUp + TreeAmp(minus)*aL_QUp) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqbg_ZH(str_) = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1
    TreeAmp_qqbg_ZH(bot_) = (TreeAmp(plus)*aR_QDn + TreeAmp(minus)*aL_QDn) * dsqrt(couplZffsq) * ZProp1


RETURN
END FUNCTION

      
      





! 
! 
! SUBROUTINE IntDip_qqb_Z_ttb(z,sHat,IDip)
! use ModParameters
! use ModKinematics
! use ModMisc
! implicit none
! real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sij,sHat
! real(8) :: dipsoft,dipfini,dipplus,epcorr
! integer :: n,emi
! 
! 
! epinv2=0d0
! epinv =0d0
! 
!    sij = sHat
! 
!    IDip(1:3) = 0d0
! do n=1,4
!       if(n.eq.1) then
!         dipsoft =ii_qq(sij,z,1) * CF
!         dipfini =ii_qq(sij,z,2) * CF
!         dipplus =ii_qq(sij,z,3) * CF
!         emi = 1
!       elseif(n.eq.2) then
!         dipsoft =ii_qq(sij,z,1) * CF
!         dipfini =ii_qq(sij,z,2) * CF
!         dipplus =ii_qq(sij,z,3) * CF
!         emi = 2
!       endif
! 
!       if(emi.eq.1) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(2) = IDip(2) + (dipfini+dipplus)
!       elseif(emi.eq.2) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(3) = IDip(3) + (dipfini+dipplus)
!       endif
!    enddo
! 
! 
! ! print *, epinv
! ! print *, "IntDip",IDip(2:3)
! 
! 
! ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
!        epcorr=epinv
! 
! !      this is for qqb-->Z-->ttb+g process
!        APsoft= 3d0/2d0*CF     * epcorr
!        APfini= (-1d0-z)*CF    * epcorr
!        APplus= 2d0*CF/(1d0-z) * epcorr
!        IDip(1) = IDip(1) + (APsoft - APplus)*2d0
!        IDip(2) = IDip(2) + (APfini + APplus)
!        IDip(3) = IDip(3) + (APfini + APplus)
! 
! ! print *, "AP",(APsoft - APplus),(APfini + APplus)
! ! print *, "sum",IDip(2:3)
! ! pause
! 
! END SUBROUTINE
! 

! 
! 
! 
! SUBROUTINE IntDip_gqb_Z_ttb(z,sHat,IDip)
! use ModParameters
! use ModKinematics
! use ModMisc
! implicit none
! real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sij,sHat
! real(8) :: dipsoft,dipfini,dipplus,epcorr
! integer :: n,emi
! 
! 
! epinv2=0d0
! epinv =0d0
! 
!    sij = sHat
! 
!    IDip(1:3) = 0d0
! 
!    dipsoft =ii_gq(sij,z,1) * TR
!    dipfini =ii_gq(sij,z,2) * TR
!    dipplus =ii_gq(sij,z,3) * TR
!    emi = 1
! 
! 
!       if(emi.eq.1) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(2) = IDip(2) + (dipfini+dipplus)
!       elseif(emi.eq.2) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(3) = IDip(3) + (dipfini+dipplus)
!       elseif(emi.eq.3) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
! !         IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
! !         IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
!       endif
! 
! 
! ! print *, epinv
! ! print *, "IntDip",IDip(2:3)
! 
! ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
!        epcorr=epinv
! 
! !      this is for gqb-->Z-->ttb+qb process
!        APsoft= 0d0
!        APfini= TR*(z**2+(1d0-z)**2) * epcorr
!        APplus= 0d0
!        IDip(1) = IDip(1) + (APsoft - APplus)
!        IDip(2) = IDip(2) + (APfini + APplus)
! 
! ! print *, "AP",(APsoft - APplus),(APfini + APplus)
! ! print *, "sum",IDip(2:3)
! ! pause
! 
! END SUBROUTINE
! 
! 
! 
! 
! SUBROUTINE IntDip_qg_Z_ttb(z,sHat,IDip)
! use ModParameters
! use ModKinematics
! use ModMisc
! implicit none
! real(8) :: IDip(1:3),APsoft,APfini,APplus,z,sij,sHat
! real(8) :: dipsoft,dipfini,dipplus,epcorr
! integer :: n,emi
! 
! 
! epinv2=0d0
! epinv =0d0
! 
!    sij = sHat
! 
!    IDip(1:3) = 0d0
! 
!    dipsoft =ii_gq(sij,z,1) * TR
!    dipfini =ii_gq(sij,z,2) * TR
!    dipplus =ii_gq(sij,z,3) * TR
!    emi = 2
! 
! 
!       if(emi.eq.1) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(2) = IDip(2) + (dipfini+dipplus)
!       elseif(emi.eq.2) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
!         IDip(3) = IDip(3) + (dipfini+dipplus)
!       elseif(emi.eq.3) then
!         IDip(1) = IDip(1) + (dipsoft-dipplus)
! !         IDip(2) = IDip(2) + (dipfini+dipplus)*0.5d0
! !         IDip(3) = IDip(3) + (dipfini+dipplus)*0.5d0
!       endif
! 
! 
! ! print *, epinv
! ! print *, "IntDip",IDip(2:3)
! 
! ! !        epcorr=epinv+2d0*dlog(renscale/facscale)
!        epcorr=epinv
! 
! !      this is for gqb-->Z-->ttb+qb process
!        APsoft= 0d0
!        APfini= TR*(z**2+(1d0-z)**2) * epcorr
!        APplus= 0d0
!        IDip(1) = IDip(1) + (APsoft - APplus)
!        IDip(3) = IDip(3) + (APfini + APplus)
! 
! ! print *, "AP",(APsoft - APplus),(APfini + APplus)
! ! print *, "sum",IDip(2:3)
! ! pause
! 
! END SUBROUTINE
! 
! 
! 




SUBROUTINE CalcFormFactor1(xe,sHat,FF1)
use ModParameters
use ModMisc
implicit none
complex(8) :: FF1
integer :: xe
real(8) :: sHat
complex(8) :: qlI2,qlI3,SI(2:3),C(1:1)


    C(1) = (-2d0,0d0)
    if( xe.ne.0 ) C(1)=(0d0,0d0)

#if useQCDLoopLib==1 || linkMELA==1
    SI(2) = qlI2(shat,0d0,0d0,Mu_Ren**2,xe)
    SI(3) = qlI3(shat,0d0,0d0,0d0,0d0,0d0,Mu_Ren**2,xe)
    FF1 = C(1) -3d0*SI(2) -2d0*sHat*SI(3)
#else
    call Error("This process requires linkQCDLoop or linkMELA enabled in the makefile")
#endif

END SUBROUTINE



SUBROUTINE Dipoles_qqb_VH(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
use ModParameters
use ModMisc
implicit none
integer :: nDipole,a,i,b,j
real(8) :: MomExt(1:4,1:7),MomExtTd(1:4,1:6),Q(1:4),QTd(1:4),KSum(1:4)
real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
real(8) :: Dipole

  if(nDipole.eq.1) then
      a=1; i=7; b=2!   initial-initial

      sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
      sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
      sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
      x = 1d0 - (sai+sbi)/sab
      v = sai/sab
      if( alpha_ii.lt.v ) then
         Dipole = 0d0
         return
      endif

      MomExtTd(1:4,a) = x*MomExt(1:4,a)
      MomExtTd(1:4,b) = MomExt(1:4,b)

      Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
      QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
      KSum(1:4) = Q(1:4)+QTd(1:4)
      
      do j=3,6
         MomExtTd(1:4,j) = MomExt(1:4,j) - 2d0*(MomExt(1:4,j).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,j).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
      enddo
      Dipole = -1d0/sai/x * 2d0*CF * (2d0/(1d0-x)-1d0-x)

  elseif(nDipole.eq.2) then
      a=2; i=7; b=1!   initial-initial

      sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
      sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
      sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
      x = 1d0 - (sai+sbi)/sab
      v = sai/sab
      if( alpha_ii.lt.v ) then
         Dipole = 0d0
         return
      endif

      MomExtTd(1:4,a) = x*MomExt(1:4,a)
      MomExtTd(1:4,b) = MomExt(1:4,b)

      Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
      QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
      KSum(1:4) = Q(1:4)+QTd(1:4)
      
      do j=3,6
        MomExtTd(1:4,j) = MomExt(1:4,j) - 2d0*(MomExt(1:4,j).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,j).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
      enddo
      Dipole = -1d0/sai/x * 2d0*CF * (2d0/(1d0-x)-1d0-x)

  endif

END SUBROUTINE




! SUBROUTINE Dipoles_gqb_Z_ttbqb(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
! use ModParameters
! use ModKinematics
! use ModMisc
! implicit none
! integer :: nDipole,a,i,b,j,k
! real(8) :: MomExt(1:4,1:11),MomExtTd(1:4,1:10),Q(1:4),QTd(1:4),KSum(1:4)
! real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
! complex(8) :: Dipole
! 
! 
!       a=1; i=5; b=2!   initial-initial
! 
!       sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
!       sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
!       sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
!       x = 1d0 - (sai+sbi)/sab
!       v = sai/sab
!       if( alpha_ii.lt.v ) then
!          Dipole = (0d0,0d0)
!          return
!       endif
! 
!       MomExtTd(1:4,a) = x*MomExt(1:4,a)
!       MomExtTd(1:4,b) = MomExt(1:4,b)
! 
!       Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
!       QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
!       KSum(1:4) = Q(1:4)+QTd(1:4)
!       MomExtTd(1:4,3) = MomExt(1:4,3) - 2d0*(MomExt(1:4,3).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,3).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
!       MomExtTd(1:4,4) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
! 
!       Dipole = -1d0/sai/x * 2d0*TR * (1d0-2d0*x*(1d0-x))
! 
! 
! END SUBROUTINE
! 
! 
! 
! 
! 
! SUBROUTINE Dipoles_qg_Z_ttbq(nDipole,MomExt,MomExtTd,Dipole)! global norm:   4d0*Pi*alpha_s
! use ModParameters
! use ModKinematics
! use ModMisc
! implicit none
! integer :: nDipole,a,i,b,j,k
! real(8) :: MomExt(1:4,1:11),MomExtTd(1:4,1:10),Q(1:4),QTd(1:4),KSum(1:4)
! real(8) :: sab,sai,sbi,sij,sik,skj,x,v,y,yp,z,Q2,mu2,mu,MomFac1,MomFac2,MomFac3
! complex(8) :: Dipole
! 
!       a=2; i=5; b=1!   initial-initial
! 
!       sab = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,b))
!       sai = 2d0*(MomExt(1:4,a).dot.MomExt(1:4,i))
!       sbi = 2d0*(MomExt(1:4,b).dot.MomExt(1:4,i))
!       x = 1d0 - (sai+sbi)/sab
!       v = sai/sab
!       if( alpha_ii.lt.v ) then
!          Dipole = (0d0,0d0)
!          return
!       endif
! 
!       MomExtTd(1:4,a) = x*MomExt(1:4,a)
!       MomExtTd(1:4,b) = MomExt(1:4,b)
! 
!       Q(1:4)   = MomExt(1:4,a)+MomExt(1:4,b)-MomExt(1:4,i)
!       QTd(1:4) = MomExtTd(1:4,a)+MomExtTd(1:4,b)
!       KSum(1:4) = Q(1:4)+QTd(1:4)
!       MomExtTd(1:4,3) = MomExt(1:4,3) - 2d0*(MomExt(1:4,3).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,3).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
!       MomExtTd(1:4,4) = MomExt(1:4,4) - 2d0*(MomExt(1:4,4).dot.KSum(1:4))/(KSum(1:4).dot.KSum(1:4))*KSum(1:4) + 2d0*(MomExt(1:4,4).dot.Q(1:4))/(Q(1:4).dot.Q(1:4))*QTd(1:4)
! 
!       Dipole = -1d0/sai/x * 2d0*TR * (1d0-2d0*x*(1d0-x))
! 
! END SUBROUTINE
! 









!     FUNCTION ii_qq(sij,x,s)
!     use ModMisc
!     use ModParameters
!     implicit none
!     real(8), intent(in)  :: x,sij
!     integer, intent(in) :: s
!     real(8) :: ii_qq
!     real(8) :: L,lx, Pqqreg
! 
! 
!         L = dlog(sij/MuRen**2)
!         lx = dlog(x)
! 
!         Pqqreg = -1d0-x
!         if (s.eq.1) then
!               ii_qq = epinv*(epinv2-L) +L**2/2d0 - Pi**2/6.0      ! CET eq.(A.4)
! !               if (scheme.eq.'fdh') ii_qq = ii_qq - 1d0/2d0
!         elseif (s.eq.2) then
!               ii_qq = -(epinv-L+lx)*Pqqreg+2d0*Pqqreg*dlog(1d0-x) - 2d0*dlog(x)/(1d0-x)+1d0-x      ! CET eq.(A.4)
!               if (alpha_ii/(1d0-x) .lt. 1d0) ii_qq=ii_qq + (2d0/(1d0-x)+Pqqreg)*dlog(alpha_ii/(1d0-x))   ! NEW
!         elseif (s.eq.3) then
!               ii_qq = 4d0*dlog(1d0-x)/(1d0-x) -(epinv-L)*2d0/(1d0-x) ! CET eq.(A.4)
!         endif
! 
! 
!     RETURN
!     END FUNCTION




! ------------------ BEGIN MCFM SM check  ------------      



      function ggHZ_pp_box(i1,i2,i3,i4,za,zb,s,mt2)
      use ModParameters
      implicit none
      complex(8)::  ggHZ_pp_box
      integer:: i1,i2,i3,i4 
      real(8):: mt2,s(1:4,1:4)
      real(8):: mZsq,mHsq,s15,s25,s12,t
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=5)
      complex(8):: qlI4,qlI3,qlI2
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5)
      complex(8):: D0(Nbox),C0(Ntri)
      complex(8):: di(Nbox),cii(Ntri)
      complex(8):: za(1:4,1:4),zb(1:4,1:4)
     

      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_pp_box=(0d0,0d0)
      
      
!======= kinematic configurations
      s12=s(i1,i2)
      mZsq=s(i3,i4) 
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4) 
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
      
    

      di(:)=(0d0,0d0)
      cii(:)=(0d0,0d0)
!
#if useQCDLoopLib==1 || linkMELA==1
     D0(d25_12)=qlI4(mZsq,0d0,0d0,mHsq,s25,s12,mt2,mt2,mt2,mt2,Mu_Ren**2,0) 
     D0(d15_12)=qlI4(mZsq,0d0,0d0,mHsq,s15,s12,mt2,mt2,mt2,mt2,Mu_Ren**2,0) 
     D0(d15_25)=qlI4(mHsq,0d0,mZsq,0d0,s15,s25,mt2,mt2,mt2,mt2,Mu_Ren**2,0) 
     C0(c25_Z)=qlI3(s25,0d0,mZsq,mt2,mt2,mt2,Mu_Ren**2,0) 
     C0(cH_25)=qlI3(mHsq,0d0,s25,mt2,mt2,mt2,Mu_Ren**2,0) 
     C0(c12)=qlI3(s12,0d0,0d0,mt2,mt2,mt2,Mu_Ren**2,0) 
     C0(c15_H)=qlI3(s15,0d0,mHsq,mt2,mt2,mt2,Mu_Ren**2,0) 
     C0(cZ_15)=qlI3(mZsq,0d0,s15,mt2,mt2,mt2,Mu_Ren**2,0) 
#else
    call Error("This process requires linkQCDLoop or linkMELA enabled in the makefile")
#endif

!------- coefficiients of integrals 

      di(d25_12)= &
      (zb(i2,i1)*(2*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
             (-2*mt2*za(i2,i3) + za(i1,i2)*za(i3,i4)*zb(i4,i1))* &
             zb(i4,i2)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) &
             - za(i1,i3)*(za(i2,i3)*zb(i3,i1) + &
               za(i2,i4)*zb(i4,i1))* &
             (t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i2) + &
               2*(2*mt2 - za(i1,i2)*zb(i2,i1))*zb(i4,i1)* &
                (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))) - &
            t(i1,i3,i4)*za(i1,i2)* &
             (za(i2,i3)*zb(i2,i1)*zb(i4,i1)* &
                (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + &
               za(i3,i4)*((za(i2,i3)*zb(i3,i1) + &
                     za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2 - &
                  zb(i4,i1)**2* &
                   (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/ &
        (2*s(i3,i4)*za(i1,i2)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))

      di(d15_12)= (zb(i2,i1)*(2*zb(i4,i1)* &
             (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
             (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))* &
             (-2*mt2*za(i1,i3) - za(i1,i2)*za(i3,i4)*zb(i4,i2)) - &
            za(i2,i3)*(za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))* &
             (t(i2,i3,i4)*za(i1,i2)*zb(i2,i1)*zb(i4,i1) + &
               2*(2*mt2 - za(i1,i2)*zb(i2,i1))* &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                zb(i4,i2)) + &
            t(i2,i3,i4)*za(i1,i2)* &
             (-(za(i1,i3)*zb(i2,i1)* &
                  (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                  zb(i4,i2)) + &
               za(i3,i4)*(-((za(i2,i3)*zb(i3,i1) + &
                       za(i2,i4)*zb(i4,i1))*zb(i4,i2)**2) + &
                  zb(i4,i1)**2* &
                   (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))))/ &
        (2*s(i3,i4)*za(i1,i2)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))


      di(d15_25)= -((4*mt2*zb(i2,i1)*(za(i1,i3)**2*zb(i3,i2)*zb(i4,i1)* &
                 (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) - &
                za(i1,i2)*za(i3,i4)* &
                 (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                 zb(i4,i2)**2 + &
                za(i1,i3)*zb(i4,i2)* &
                 (za(i2,i3)**2*zb(i3,i1)*zb(i3,i2) + &
                   za(i2,i4)*zb(i4,i1)* &
                    (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + &
                   za(i2,i3)* &
                    (za(i1,i4)*zb(i3,i1)*zb(i4,i1) + &
                      za(i2,i4)* &
                       (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) &
                      ))))/ &
            (za(i1,i2)*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
              (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))) + &
           (2*za(i1,i3)*(za(i1,i3)**2*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                  (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                  (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2)) + &
                 za(i1,i3)*(za(i2,i3)**3*zb(i3,i1)*zb(i3,i2)**2* &
                     (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) + &
                    za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)* &
                     (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2))* &
                     (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) &
                     + za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)* &
                     (za(i3,i4)*zb(i3,i2)*zb(i4,i3) + &
                       za(i1,i4)* &
                        (2*zb(i3,i2)*zb(i4,i1) + &
                          zb(i2,i1)*zb(i4,i3))) + &
                    za(i2,i3)*za(i2,i4)* &
                     (za(i2,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                        (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) &
                        + za(i1,i4)*zb(i4,i1)* &
                        (2*zb(i3,i2)*zb(i4,i1) + &
                           zb(i2,i1)*zb(i4,i3))**2 - &
                       za(i3,i4)*zb(i4,i3)* &
                        (-(zb(i3,i2)**2*zb(i4,i1)**2) + &
                          zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + &
                          zb(i2,i1)**2*zb(i4,i3)**2))) + &
                 za(i1,i4)*(za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)**2* &
                     (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + &
                    za(i2,i3)**3*zb(i3,i2)* &
                     (2*zb(i3,i2)**2*zb(i4,i1)**2 + &
                       3*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + &
                       zb(i2,i1)**2*zb(i4,i3)**2) + &
                    za(i2,i3)*za(i2,i4)*zb(i4,i1)*zb(i4,i2)* &
                     (-(za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) + &
                       za(i1,i4)*zb(i4,i1)* &
                        (2*zb(i3,i2)*zb(i4,i1) + &
                          zb(i2,i1)*zb(i4,i3)) + &
                       2*za(i2,i4)*zb(i4,i2)* &
                        (2*zb(i3,i2)*zb(i4,i1) + &
                          zb(i2,i1)*zb(i4,i3))) + &
                    za(i2,i3)**2* &
                     (za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2* &
                        (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) &
                        + za(i3,i4)*zb(i4,i3)* &
                        (-(zb(i3,i2)**2*zb(i4,i1)**2) + &
                          zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3) + &
                          zb(i2,i1)**2*zb(i4,i3)**2) + &
                       za(i2,i4)*zb(i4,i2)* &
                        (5*zb(i3,i2)**2*zb(i4,i1)**2 + &
                          5*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)* &
                           zb(i4,i3) + zb(i2,i1)**2*zb(i4,i3)**2)))) &
                + za(i1,i2)* &
               (za(i1,i4)*zb(i4,i2)* &
                  (-(za(i2,i4)*za(i3,i4)*zb(i4,i1)**2*zb(i4,i2)* &
                       (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2))) &
                     + za(i2,i3)*zb(i4,i1)*zb(i4,i2)* &
                     (za(i1,i4)* &
                        (za(i2,i4)*zb(i2,i1) - za(i3,i4)*zb(i3,i1))* &
                        zb(i4,i1) - &
                       2*za(i2,i4)*za(i3,i4)* &
                        (zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3))) &
                      - za(i2,i3)**2* &
                     (-(za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)* &
                          zb(i4,i2)) + &
                       za(i3,i4)* &
                        (zb(i3,i2)*zb(i4,i1) + &
                           zb(i2,i1)*zb(i4,i3))**2)) + &
                 za(i1,i3)**2*zb(i3,i2)* &
                  (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                  (za(i2,i3)*zb(i2,i1)* &
                     (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)) &
                     + zb(i4,i1)* &
                     (za(i2,i4)*zb(i2,i1)*zb(i4,i2) - &
                       za(i3,i4)* &
                        (zb(i3,i2)*zb(i4,i1) + &
                          2*zb(i2,i1)*zb(i4,i3)))) + &
                 za(i1,i3)*(za(i2,i4)*zb(i4,i1)**2*zb(i4,i2)* &
                     (za(i1,i4)* &
                        (za(i2,i4)*zb(i2,i1) - &
                          2*za(i3,i4)*zb(i3,i1))*zb(i4,i2) + &
                       za(i3,i4)*zb(i3,i2)* &
                        (za(i2,i4)*zb(i4,i2) - &
                          2*za(i3,i4)*zb(i4,i3))) + &
                    za(i2,i3)**2* &
                     (za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i2)* &
                        (3*zb(i3,i2)*zb(i4,i1) + &
                          zb(i2,i1)*zb(i4,i3)) + &
                       za(i3,i4)* &
                        (zb(i3,i2)**3*zb(i4,i1)**2 - &
                          zb(i2,i1)**2*zb(i3,i2)*zb(i4,i3)**2)) + &
                    za(i2,i3)* &
                     (za(i1,i4)*zb(i4,i1)*zb(i4,i2)* &
                        (-2*za(i3,i4)*zb(i3,i1)**2*zb(i4,i2) + &
                          za(i2,i4)*zb(i2,i1)* &
                           (3*zb(i3,i2)*zb(i4,i1) + &
                             zb(i3,i1)*zb(i4,i2) + &
                             zb(i2,i1)*zb(i4,i3))) + &
                       2*za(i3,i4)* &
                        (za(i2,i4)*zb(i3,i2)**2*zb(i4,i1)**2* &
                           zb(i4,i2) + &
                          za(i3,i4)*zb(i4,i3)* &
                           (-(zb(i3,i2)**2*zb(i4,i1)**2) + &
                             zb(i2,i1)*zb(i3,i2)*zb(i4,i1)* &
                             zb(i4,i3) + zb(i2,i1)**2*zb(i4,i3)**2)) &
                       ))))/ &
            (za(i1,i2)**2*(za(i2,i3)*zb(i3,i1) + &
                za(i2,i4)*zb(i4,i1))* &
              (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))/ &
        (2*s(i3,i4))


!---- triangle coefficiients
      cii(c25_Z)=(za(i1,i3)*(za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))* &
          (za(i1,i3)*zb(i3,i2)*zb(i4,i1)* &
             (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + &
            zb(i4,i2)*(za(i2,i3)**2*zb(i3,i1)*zb(i3,i2) + &
               za(i1,i2)*zb(i2,i1)* &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + &
               za(i2,i4)*zb(i4,i1)* &
                (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + &
               za(i2,i3)*(za(i1,i4)*zb(i3,i1)*zb(i4,i1) + &
                  za(i2,i4)* &
                   (2*zb(i3,i2)*zb(i4,i1) + zb(i2,i1)*zb(i4,i3)))))) &
         /(za(i1,i2)**2*za(i3,i4)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      cii(cZ_15)= (za(i2,i3)*(za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))* &
          (za(i2,i3)*zb(i3,i1)*zb(i4,i2)* &
             (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + &
            zb(i4,i1)*(za(i1,i3)**2*zb(i3,i1)*zb(i3,i2) + &
               za(i1,i2)*zb(i2,i1)* &
                (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + &
               za(i1,i4)*zb(i4,i2)* &
                (za(i1,i4)*zb(i4,i1) + za(i2,i4)*zb(i4,i2)) + &
               za(i1,i3)*(za(i2,i4)*zb(i3,i2)*zb(i4,i2) + &
                  za(i1,i4)* &
                   (2*zb(i3,i1)*zb(i4,i2) - zb(i2,i1)*zb(i4,i3)))))) &
         /(za(i1,i2)**2*za(i3,i4)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
          (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))*zb(i4,i3))

      cii(cH_25)=-((za(i2,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i1) + &
              (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))* &
               (za(i1,i3)*zb(i3,i1)*zb(i4,i1) + &
                 za(i1,i4)*zb(i4,i1)**2 + &
                 (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                  zb(i4,i2)) + &
              za(i1,i2)*zb(i2,i1)* &
               (za(i1,i3)*zb(i3,i1)*zb(i4,i1) + &
                 za(i1,i4)*zb(i4,i1)**2 + &
                 2*za(i2,i3)*zb(i3,i1)*zb(i4,i2) + &
                 2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) - &
                 za(i2,i3)*zb(i2,i1)*zb(i4,i3))))/ &
          (s(i3,i4)*za(i1,i2)**2* &
            (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))))

      cii(c15_H)= -((za(i1,i3)*(za(i1,i2)**2*zb(i2,i1)**2*zb(i4,i2) + &
              (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))* &
               (za(i2,i3)*zb(i3,i2)*zb(i4,i2) + &
                 za(i2,i4)*zb(i4,i2)**2 + &
                 zb(i4,i1)*(za(i1,i3)*zb(i3,i2) + &
                    za(i1,i4)*zb(i4,i2))) + &
              za(i1,i2)*zb(i2,i1)* &
               (2*za(i1,i3)*zb(i3,i2)*zb(i4,i1) + &
                 za(i2,i3)*zb(i3,i2)*zb(i4,i2) + &
                 2*za(i1,i4)*zb(i4,i1)*zb(i4,i2) + &
                 za(i2,i4)*zb(i4,i2)**2 + &
                 za(i1,i3)*zb(i2,i1)*zb(i4,i3))))/ &
          (s(i3,i4)*za(i1,i2)**2* &
            (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))

      cii(c12)=  -((zb(i2,i1)*((zb(i4,i1)* &
                 (-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1)))/ &
               (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) - &
              (zb(i4,i2)*(za(i1,i3)*zb(i2,i1) + &
                   za(i3,i4)*zb(i4,i2)))/ &
               (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2))))/ &
          s(i3,i4))



      ggHZ_pp_box=(0d0,0d0)
        do i=1,Nbox
           ggHZ_pp_box=ggHZ_pp_box+D0(i)*di(i)
        enddo
        do i=1,Ntri
           ggHZ_pp_box=ggHZ_pp_box+C0(i)*cii(i)
        enddo

      return
      end



      
      function ggHZ_mp_box(i1,i2,i3,i4,za,zb,s,mt2)
      implicit none
      complex(8)::  ggHZ_mp_box,za(1:4,1:4),zb(1:4,1:4)
      integer:: i1,i2,i3,i4
      real(8):: mt2,s(1:4,1:4)
      real(8):: mZsq,mHsq,s15,s25,s12,t
      integer:: Nbox,Ntri,i
      parameter(Nbox=3,Ntri=6)

      complex(8):: qlI4,qlI3,qlI2
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15,c12_Z_H
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5,c12_Z_H=6)
!--- nb I swapped the notation wrt to KC for 3m triangle, so that I
!---- can unify two rotuines and save calls to QCDLoop.
      complex(8):: D0(Nbox),C0(Ntri)
      complex(8):: di(Nbox),cii(Ntri)

!       complex(8)::  ggHZ_mp_2me_b1,ggHZ_mp_2mh_b1
!       complex(8)::  ggHZ_mp_tri1,ggHZ_mp_tri2,ggHZ_mp_3mtri


      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_box=(0d0,0d0)


!======= kinematic configurations
      s12=s(i1,i2)
      mZsq=s(i3,i4)
      s25=t(i1,i3,i4)
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4)
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

      di(:)=(0d0,0d0)
      cii(:)=(0d0,0d0)


!----- fill integrals from QCD loop : boxes

#if useQCDLoopLib==1 || linkMELA==1
     D0(d25_12)=qlI4(mZsq,0d0,0d0,mHsq,s25,s12,mt2,mt2,mt2,mt2,Mu_Ren**2,0)
     D0(d15_12)=qlI4(mZsq,0d0,0d0,mHsq,s15,s12,mt2,mt2,mt2,mt2,Mu_Ren**2,0)
     D0(d15_25)=qlI4(mHsq,0d0,mZsq,0d0,s15,s25,mt2,mt2,mt2,mt2,Mu_Ren**2,0)

! ----- fill integrals from QCD loop : triangles
     C0(c25_Z)=qlI3(s25,0d0,mZsq,mt2,mt2,mt2,Mu_Ren**2,0)
     C0(cH_25)=qlI3(mHsq,0d0,s25,mt2,mt2,mt2,Mu_Ren**2,0)
     C0(c12)=qlI3(s12,0d0,0d0,mt2,mt2,mt2,Mu_Ren**2,0)
     C0(c15_H)=qlI3(s15,0d0,mHsq,mt2,mt2,mt2,Mu_Ren**2,0)
     C0(cZ_15)=qlI3(mZsq,0d0,s15,mt2,mt2,mt2,Mu_Ren**2,0)
     C0(c12_Z_H)=qlI3(s12,mZsq,mHsq,mt2,mt2,mt2,Mu_Ren**2,0)
#else
    call Error("This process requires linkQCDLoop or linkMELA enabled in the makefile")
#endif

!------ box coefficiients
      di(d25_12)= ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,s,mt2)
      di(d15_12)= -ggHZ_mp_2mh_b1(i2,i1,i4,i3,zb,za,s,mt2)
      di(d15_25)= ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,s,mt2)


!==== triangle coefficiients

      cii(c25_Z)=ggHZ_mp_tri2(i1,i2,i3,i4,za,zb,s,mt2)
      cii(cZ_15)=-ggHZ_mp_tri2(i2,i1,i4,i3,zb,za,s,mt2)
      cii(cH_25)=-ggHZ_mp_tri1(i1,i2,i3,i4,za,zb,s,mt2)
      cii(c15_H)=ggHZ_mp_tri1(i2,i1,i4,i3,zb,za,s,mt2)
      cii(c12_Z_H)=ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb,s)
      cii(c12)= (-(za(i1,i3)*(za(i1,i4)*za(i2,i3) + &
               za(i1,i3)*za(i2,i4))*zb(i2,i1)* &
             (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
             zb(i4,i3)) + &
          za(i1,i2)**2*za(i3,i4)*zb(i2,i1)* &
           (-(za(i2,i3)*zb(i2,i1)) + za(i3,i4)*zb(i4,i1))* &
           zb(i4,i3) + za(i1,i2)* &
           (2*za(i1,i3)**2*za(i2,i4)*zb(i2,i1)*zb(i3,i1)* &
              zb(i4,i1) + &
             2*za(i2,i4)**2*za(i3,i4)*zb(i4,i1)* &
              zb(i4,i2)**2 - &
             za(i2,i3)**2*zb(i2,i1)* &
              (2*za(i2,i4)*zb(i3,i2)*zb(i4,i2) + &
                za(i1,i4)*zb(i2,i1)*zb(i4,i3)) + &
             za(i1,i3)*zb(i2,i1)* &
              (2*za(i1,i4)*za(i2,i4)*zb(i4,i1)**2 + &
                2*za(i2,i4)**2*zb(i4,i1)*zb(i4,i2) - &
                za(i2,i3)*za(i3,i4)*zb(i3,i1)*zb(i4,i3) + &
                za(i2,i4)* &
                 (-(za(i2,i3)*zb(i2,i1)) + &
                   2*za(i3,i4)*zb(i4,i1))*zb(i4,i3)) + &
             za(i2,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                 (-2*za(i2,i4)*zb(i4,i2) + &
                   za(i3,i4)*zb(i4,i3)) - &
                2*za(i2,i4)*zb(i4,i2)* &
                 (za(i2,i4)*zb(i2,i1)*zb(i4,i2) + &
                   za(i3,i4)* &
                    (-(zb(i3,i1)*zb(i4,i2)) + &
                      2*zb(i2,i1)*zb(i4,i3))))))/ &
        (2.*za(i2,i4)*za(i3,i4)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2* &
          zb(i4,i3))


      ggHZ_mp_box=(0d0,0d0)
!===== make amplitude
        do i=1,Nbox
           ggHZ_mp_box=ggHZ_mp_box+D0(i)*di(i)
        enddo
        do i=1,Ntri
           ggHZ_mp_box=ggHZ_mp_box+C0(i)*cii(i)
        enddo

      return
      end


      function ggHZ_mp_2mh_b1(i1,i2,i3,i4,za,zb,s,mt2)
      implicit none
      complex(8):: ggHZ_mp_2mh_b1,za(1:4,1:4),zb(1:4,1:4)
!---- coefficiient of 2mh box
      integer:: i1,i2,i3,i4
      real(8):: mt2,t,s(1:4,1:4)
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

        ggHZ_mp_2mh_b1=(-2*za(i1,i3)*(2*mt2*zb(i4,i1)* &
              (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
              (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)) + &
             t(i1,i3,i4)*za(i1,i2)*zb(i2,i1)* &
              (t(i1,i3,i4)*zb(i4,i1) + &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i2))) + &
          zb(i4,i2)*(-(t(i1,i3,i4)*za(i1,i2)*za(i3,i4)* &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i2)) + &
             2*za(i2,i3)*(t(i1,i3,i4)**2*za(i1,i2)*zb(i2,i1) + &
                2*mt2*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                 (za(i1,i3)*zb(i3,i2) + za(i1,i4)*zb(i4,i2)))) + &
          t(i1,i3,i4)*za(i1,i3)**2*zb(i2,i1)* &
           (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))*zb(i4,i3))/ &
        (2.*s(i3,i4)*(za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2)

        return
        end

      function ggHZ_mp_2me_b1(i1,i2,i3,i4,za,zb,s,mt2)
      implicit none
      complex(8):: ggHZ_mp_2me_b1,za(1:4,1:4),zb(1:4,1:4)
!---- coefficiient of 2mh box
      integer:: i1,i2,i3,i4
      real(8):: mt2,s(1:4,1:4)

      ggHZ_mp_2me_b1= &
      (-(za(i1,i2)**3*za(i3,i4)**2*zb(i2,i1)**2*zb(i4,i3)* &
             (za(i2,i3)**2*zb(i3,i1)**2*zb(i3,i2)*zb(i4,i1)* &
                zb(i4,i2) + &
               2*za(i2,i3)*zb(i3,i1)*zb(i4,i1)* &
                (za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2 + &
                  zb(i4,i3)* &
                   (za(i1,i3)*zb(i2,i1)*zb(i3,i1) + &
                     za(i3,i4)* &
                      (-(zb(i3,i1)*zb(i4,i2)) + &
                        zb(i2,i1)*zb(i4,i3)))) + &
               za(i2,i4)*zb(i4,i1)**2* &
                (za(i2,i4)*zb(i4,i2)* &
                   (zb(i3,i1)*zb(i4,i2) + zb(i2,i1)*zb(i4,i3)) &
                    - 2*zb(i4,i3)* &
                   (-(za(i1,i3)*zb(i2,i1)*zb(i3,i1)) + &
                     za(i3,i4)* &
                      (zb(i3,i1)*zb(i4,i2) + &
                        zb(i2,i1)*zb(i4,i3)))))) - &
          (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1))* &
           (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))* &
           (-2*za(i1,i4)**2*za(i2,i3)*zb(i4,i1)**2* &
              (za(i2,i3)**2*zb(i3,i1)*zb(i3,i2)**2* &
                 zb(i4,i1) + &
                za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)* &
                 zb(i4,i2)**2 - &
                za(i2,i3)*za(i2,i4)* &
                 (zb(i3,i2)**2*zb(i4,i1)**2 - &
                   4*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i2) + zb(i3,i1)**2*zb(i4,i2)**2)) + &
             za(i1,i3)**2*za(i2,i4)*zb(i3,i1)* &
              (-2*za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)* &
                 (zb(i3,i2)*zb(i4,i1) - 2*zb(i3,i1)*zb(i4,i2)) &
                  + za(i2,i3)**2*zb(i3,i1)* &
                 (zb(i3,i2)**2*zb(i4,i1)**2 + &
                   zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + &
                   zb(i2,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) + &
                za(i2,i3)*za(i2,i4)*zb(i4,i1)* &
                 (-3*zb(i3,i2)**2*zb(i4,i1)**2 + &
                   zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + &
                   zb(i3,i2)*zb(i4,i1)* &
                    (6*zb(i3,i1)*zb(i4,i2) + &
                      zb(i2,i1)*zb(i4,i3)))) + &
             za(i1,i3)*za(i1,i4)*zb(i4,i1)* &
              (-2*za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)* &
                 zb(i4,i2) + &
                2*za(i2,i4)**3*zb(i3,i1)*zb(i4,i1)**2* &
                 zb(i4,i2)**2 + &
                za(i2,i3)**2*za(i2,i4)*zb(i3,i1)* &
                 (7*zb(i3,i2)**2*zb(i4,i1)**2 + &
                   zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2) + &
                   zb(i3,i2)*zb(i4,i1)* &
                    (-10*zb(i3,i1)*zb(i4,i2) + &
                      zb(i2,i1)*zb(i4,i3))) - &
                za(i2,i3)*za(i2,i4)**2*zb(i4,i1)* &
                 (zb(i3,i2)**2*zb(i4,i1)**2 + &
                   zb(i3,i1)*zb(i4,i2)* &
                    (5*zb(i3,i1)*zb(i4,i2) + &
                      zb(i2,i1)*zb(i4,i3)) - &
                   zb(i3,i2)*zb(i4,i1)* &
                    (8*zb(i3,i1)*zb(i4,i2) + &
                      zb(i2,i1)*zb(i4,i3))))) + &
          za(i1,i2)*(za(i1,i3)**3*zb(i2,i1)*zb(i3,i1)**2* &
              (za(i2,i3)**2*zb(i3,i1)**2 - &
                za(i2,i4)**2*zb(i4,i1)**2)* &
              (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))* &
              (-(zb(i3,i2)*zb(i4,i1)) + zb(i3,i1)*zb(i4,i2)) &
              - za(i1,i3)**2* &
              (za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)* &
                 (3*za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                    (zb(i3,i2)*zb(i4,i1) - &
                      zb(i3,i1)*zb(i4,i2)) + &
                   zb(i4,i2)* &
                    (-2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)* &
                       zb(i4,i1) + &
                      za(i2,i4)*zb(i2,i1)* &
                       (zb(i3,i2)*zb(i4,i1) + &
                         zb(i3,i1)*zb(i4,i2)))) + &
                za(i2,i3)*za(i2,i4)**2*zb(i3,i1)*zb(i4,i1)* &
                 (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                    (zb(i3,i2)**2*zb(i4,i1)**2 - &
                      3*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i2) - 4*zb(i3,i1)**2*zb(i4,i2)**2 &
                      ) - &
                   2*zb(i2,i1)**2*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i3)*(-2*mt2 + za(i3,i4)*zb(i4,i3)) &
                    + zb(i2,i1)* &
                    (zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (-8*mt2 + 3*za(i1,i4)*zb(i4,i1) + &
                         7*za(i2,i4)*zb(i4,i2) - &
                         17*za(i3,i4)*zb(i4,i3)) - &
                      zb(i3,i2)**2*zb(i4,i1)**2* &
                       (-8*mt2 + za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2) - &
                         4*za(i3,i4)*zb(i4,i3)) + &
                      2*zb(i3,i1)**2*zb(i4,i2)**2* &
                       (4*mt2 - za(i1,i4)*zb(i4,i1) + &
                         2*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i4)**3*zb(i4,i1)**2* &
                 (za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2* &
                    (zb(i3,i2)*zb(i4,i1) - &
                      3*zb(i3,i1)*zb(i4,i2)) + &
                   zb(i2,i1)**2*zb(i4,i3)* &
                    (2*mt2*zb(i3,i2)*zb(i4,i1) - &
                      za(i3,i4)*zb(i3,i1)*zb(i4,i2)*zb(i4,i3)) &
                     + zb(i2,i1)* &
                    (-2*mt2*zb(i3,i2)**2*zb(i4,i1)**2 + &
                      zb(i3,i1)**2*zb(i4,i2)**2* &
                       (4*mt2 + za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2) - &
                         8*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (2*mt2 - za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2) + &
                         3*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)**2*za(i2,i4)*zb(i3,i1)**2* &
                 (-(za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2* &
                      (5*zb(i3,i2)*zb(i4,i1) + &
                        zb(i3,i1)*zb(i4,i2))) + &
                   zb(i2,i1)**2*zb(i4,i3)* &
                    (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                       zb(i4,i3) + &
                      2*zb(i3,i2)*zb(i4,i1)* &
                       (mt2 - za(i3,i4)*zb(i4,i3))) + &
                   zb(i2,i1)* &
                    (-(zb(i3,i1)**2*zb(i4,i2)**2* &
                         (-4*mt2 + 3*za(i1,i4)*zb(i4,i1) + &
                           za(i2,i4)*zb(i4,i2))) + &
                      2*zb(i3,i2)**2*zb(i4,i1)**2* &
                       (5*mt2 + za(i1,i4)*zb(i4,i1) - &
                         4*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (-10*mt2 + za(i1,i4)*zb(i4,i1) + &
                         7*za(i2,i4)*zb(i4,i2) + &
                         5*za(i3,i4)*zb(i4,i3))))) + &
             za(i1,i4)**2*zb(i4,i1)* &
              (2*za(i2,i4)**3*za(i3,i4)*zb(i3,i1)* &
                 zb(i4,i1)**3*zb(i4,i2)**3 + &
                za(i2,i3)**4*zb(i2,i1)*zb(i3,i1)**2*zb(i3,i2)* &
                 zb(i4,i2)* &
                 (zb(i3,i2)*zb(i4,i1) + zb(i3,i1)*zb(i4,i2)) &
                 + za(i2,i3)**3*zb(i3,i1)* &
                 (2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)**2* &
                    zb(i4,i1)**2*zb(i4,i2) - &
                   zb(i2,i1)* &
                    (zb(i3,i1)**2*zb(i4,i2)**2* &
                       (4*mt2 + za(i2,i4)*zb(i4,i2)) - &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (8*mt2 + 7*za(i2,i4)*zb(i4,i2)) + &
                      6*za(i3,i4)*zb(i3,i2)**2*zb(i4,i1)**2* &
                       zb(i4,i3))) + &
                za(i2,i3)*za(i2,i4)**2*zb(i4,i1)**2*zb(i4,i2)* &
                 (2*za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                    (2*zb(i3,i2)*zb(i4,i1) + &
                      zb(i3,i1)*zb(i4,i2)) + &
                   zb(i2,i1)* &
                    (zb(i3,i1)*zb(i4,i2)* &
                       (8*mt2 + za(i2,i4)*zb(i4,i2) - &
                         10*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i2)*zb(i4,i1)* &
                       (-4*mt2 + za(i2,i4)*zb(i4,i2) + &
                         4*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)**2*za(i2,i4)*zb(i4,i1)* &
                 (2*za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i2)* &
                    (zb(i3,i2)*zb(i4,i1) + &
                      2*zb(i3,i1)*zb(i4,i2)) + &
                   zb(i2,i1)* &
                    (zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (4*mt2 + 7*za(i2,i4)*zb(i4,i2) - &
                         20*za(i3,i4)*zb(i4,i3)) + &
                      2*zb(i3,i1)**2*zb(i4,i2)**2* &
                       (2*mt2 + za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i2)**2*zb(i4,i1)**2* &
                       (-(za(i2,i4)*zb(i4,i2)) + &
                         6*za(i3,i4)*zb(i4,i3))))) + &
             za(i1,i3)*za(i1,i4)* &
              (za(i2,i3)**4*zb(i2,i1)*zb(i3,i1)**3*zb(i3,i2)* &
                 zb(i4,i2)* &
                 (zb(i3,i2)*zb(i4,i1) + zb(i3,i1)*zb(i4,i2)) &
                 - za(i2,i4)**3*zb(i4,i1)**3*zb(i4,i2)* &
                 (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                    (zb(i3,i2)*zb(i4,i1) - &
                      5*zb(i3,i1)*zb(i4,i2)) - &
                   za(i3,i4)*zb(i2,i1)**2*zb(i4,i3)**2 + &
                   zb(i2,i1)* &
                    (zb(i3,i1)*zb(i4,i2)* &
                       (8*mt2 + za(i2,i4)*zb(i4,i2) - &
                         6*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i2)*zb(i4,i1)* &
                       (-4*mt2 + za(i2,i4)*zb(i4,i2) + &
                         za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)**3*zb(i3,i1)**2* &
                 (za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i2)* &
                    (3*zb(i3,i2)*zb(i4,i1) + &
                      zb(i3,i1)*zb(i4,i2)) + &
                   2*mt2*zb(i2,i1)**2*zb(i3,i1)*zb(i4,i2)* &
                    zb(i4,i3) + &
                   zb(i2,i1)* &
                    (zb(i3,i1)**2*zb(i4,i2)**2* &
                       (2*mt2 - za(i2,i4)*zb(i4,i2)) + &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (-6*mt2 + 2*za(i1,i4)*zb(i4,i1) + &
                         6*za(i2,i4)*zb(i4,i2) - &
                         3*za(i3,i4)*zb(i4,i3)) - &
                      zb(i3,i2)**2*zb(i4,i1)**2* &
                       (-8*mt2 + 2*za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2) + &
                         4*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)**2*za(i2,i4)*zb(i3,i1)*zb(i4,i1)* &
                 (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                    (2*zb(i3,i2)**2*zb(i4,i1)**2 + &
                      9*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i2) + zb(i3,i1)**2*zb(i4,i2)**2) &
                    + zb(i2,i1)**2*zb(i4,i3)* &
                    (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i3) + &
                      zb(i3,i1)*zb(i4,i2)* &
                       (4*mt2 - za(i3,i4)*zb(i4,i3))) + &
                   zb(i2,i1)* &
                    (zb(i3,i1)**2*zb(i4,i2)**2* &
                       (8*mt2 + 2*za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2)) + &
                      zb(i3,i2)**2*zb(i4,i1)**2* &
                       (4*mt2 - 2*za(i1,i4)*zb(i4,i1) - &
                         za(i2,i4)*zb(i4,i2) + &
                         18*za(i3,i4)*zb(i4,i3)) - &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (8*mt2 + 27*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)*za(i2,i4)**2*zb(i4,i1)**2* &
                 (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                    (-(zb(i3,i2)**2*zb(i4,i1)**2) + &
                      7*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i2) + 6*zb(i3,i1)**2*zb(i4,i2)**2 &
                      ) + &
                   2*zb(i2,i1)**2*zb(i4,i3)* &
                    (mt2*zb(i3,i1)*zb(i4,i2) + &
                      za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) &
                     + zb(i2,i1)* &
                    (zb(i3,i2)**2*zb(i4,i1)**2* &
                       (-4*mt2 + za(i2,i4)*zb(i4,i2) - &
                         2*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i3,i2)*zb(i4,i1)*zb(i4,i2)* &
                       (2*mt2 - 2*za(i1,i4)*zb(i4,i1) - &
                         6*za(i2,i4)*zb(i4,i2) + &
                         23*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)**2*zb(i4,i2)**2* &
                       (2*za(i1,i4)*zb(i4,i1) + &
                         za(i2,i4)*zb(i4,i2) - &
                         2*(mt2 + 9*za(i3,i4)*zb(i4,i3))))))) &
           + za(i1,i2)**2*zb(i2,i1)* &
           (-(za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)* &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))* &
                (za(i2,i3)*zb(i3,i1)* &
                   (za(i3,i4)*zb(i3,i1)*zb(i4,i2)*zb(i4,i3) + &
                     zb(i3,i2)*zb(i4,i1)* &
                      (4*mt2 - 3*za(i3,i4)*zb(i4,i3))) + &
                  za(i2,i4)*zb(i4,i1)* &
                   (-3*za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                      zb(i4,i3) + &
                     zb(i3,i2)*zb(i4,i1)* &
                      (4*mt2 + za(i3,i4)*zb(i4,i3))))) + &
             za(i1,i4)*(za(i2,i4)**2*za(i3,i4)*zb(i4,i1)**3* &
                 zb(i4,i2)* &
                 (zb(i3,i1)*zb(i4,i2)* &
                    (4*mt2 + za(i2,i4)*zb(i4,i2) - &
                      4*za(i3,i4)*zb(i4,i3)) - &
                   zb(i2,i1)*zb(i4,i3)* &
                    (-4*mt2 + za(i2,i4)*zb(i4,i2) + &
                      2*za(i3,i4)*zb(i4,i3))) + &
                za(i2,i3)**2*zb(i3,i1)*zb(i4,i1)* &
                 (2*za(i3,i4)* &
                    (2*mt2*zb(i3,i1)**2*zb(i4,i2)**2 - &
                      2*zb(i3,i1)* &
                       (2*mt2*zb(i2,i1) + &
                         za(i3,i4)*zb(i3,i2)*zb(i4,i1))* &
                       zb(i4,i2)*zb(i4,i3) + &
                      3*za(i3,i4)*zb(i2,i1)*zb(i3,i2)* &
                       zb(i4,i1)*zb(i4,i3)**2) + &
                   za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2* &
                    (za(i3,i4)* &
                       (2*zb(i3,i2)*zb(i4,i1) + &
                         zb(i3,i1)*zb(i4,i2)) + &
                      zb(i2,i1)* &
                       (8*mt2 - 7*za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)**3*zb(i3,i1)**2*zb(i4,i2)* &
                 (za(i3,i4)*zb(i3,i1)*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i2) + &
                   zb(i2,i1)* &
                    (-2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i3) + &
                      zb(i3,i1)*zb(i4,i2)* &
                       (4*mt2 - za(i3,i4)*zb(i4,i3)))) + &
                za(i2,i3)*za(i2,i4)*zb(i4,i1)**2* &
                 (za(i2,i4)*zb(i4,i2)* &
                    (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                       (zb(i3,i2)*zb(i4,i1) + &
                         2*zb(i3,i1)*zb(i4,i2)) + &
                      zb(i2,i1)* &
                       (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)* &
                          zb(i4,i3) + &
                         zb(i3,i1)*zb(i4,i2)* &
                          (4*mt2 - 7*za(i3,i4)*zb(i4,i3)))) - &
                   2*za(i3,i4)* &
                    (3*za(i3,i4)*zb(i2,i1)*zb(i3,i2)* &
                       zb(i4,i1)*zb(i4,i3)**2 + &
                      2*zb(i3,i1)**2*zb(i4,i2)**2* &
                       (-2*mt2 + za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i4,i2)*zb(i4,i3)* &
                       (2*za(i3,i4)*zb(i3,i2)*zb(i4,i1) + &
                         zb(i2,i1)* &
                          (2*mt2 - 5*za(i3,i4)*zb(i4,i3)))))) &
              + za(i1,i3)* &
              (za(i2,i3)**3*zb(i3,i1)**3*zb(i3,i2)*zb(i4,i2)* &
                 (za(i3,i4)*zb(i3,i1)*zb(i4,i2) + &
                   zb(i2,i1)*(4*mt2 - za(i3,i4)*zb(i4,i3))) + &
                za(i2,i3)**2*zb(i3,i1)**2* &
                 (-2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                    (-2*za(i3,i4)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i3) + &
                      zb(i3,i1)*zb(i4,i2)* &
                       (2*mt2 + za(i3,i4)*zb(i4,i3))) + &
                   za(i3,i4)* &
                    (-(za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2* &
                         zb(i4,i3)) + &
                      2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)* &
                       zb(i4,i3)* &
                       (-4*mt2 + za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i4,i2)* &
                       (zb(i3,i2)*zb(i4,i1)* &
                          (4*mt2 - 3*za(i3,i4)*zb(i4,i3)) + &
                         zb(i2,i1)*zb(i4,i3)* &
                          (4*mt2 + za(i3,i4)*zb(i4,i3)))) + &
                   za(i2,i4)*zb(i4,i2)* &
                    (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                       (2*zb(i3,i2)*zb(i4,i1) + &
                         zb(i3,i1)*zb(i4,i2)) + &
                      zb(i2,i1)* &
                       (-2*za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                          zb(i4,i3) + &
                         zb(i3,i2)*zb(i4,i1)* &
                          (8*mt2 + za(i3,i4)*zb(i4,i3))))) + &
                za(i2,i4)**2*zb(i4,i1)**2* &
                 (2*za(i1,i4)*zb(i2,i1)*zb(i3,i1)*zb(i4,i1)* &
                    zb(i4,i2)*(-2*mt2 + za(i3,i4)*zb(i4,i3)) &
                    + za(i3,i4)* &
                    (zb(i3,i1)**2*zb(i4,i2)**2* &
                       (za(i2,i4)*zb(i4,i2) - &
                         5*za(i3,i4)*zb(i4,i3)) + &
                      zb(i3,i1)*zb(i4,i2)* &
                       (2*zb(i2,i1)*zb(i4,i3)* &
                          (2*mt2 + za(i2,i4)*zb(i4,i2) - &
                           3*za(i3,i4)*zb(i4,i3)) + &
                         zb(i3,i2)*zb(i4,i1)* &
                          (4*mt2 + za(i3,i4)*zb(i4,i3))) + &
                      zb(i2,i1)*zb(i4,i3)* &
                       (-(za(i3,i4)*zb(i2,i1)*zb(i4,i3)**2) + &
                         zb(i3,i2)*zb(i4,i1)* &
                          (4*mt2 - za(i2,i4)*zb(i4,i2) + &
                           za(i3,i4)*zb(i4,i3))))) + &
                za(i2,i3)*za(i2,i4)*zb(i3,i1)*zb(i4,i1)* &
                 (4*za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                    (-2*mt2*zb(i3,i1)*zb(i4,i2) + &
                      za(i3,i4)*zb(i3,i2)*zb(i4,i1)*zb(i4,i3)) &
                     + za(i2,i4)*zb(i4,i2)* &
                    (za(i3,i4)*zb(i3,i1)*zb(i4,i2)* &
                       (zb(i3,i2)*zb(i4,i1) + &
                         2*zb(i3,i1)*zb(i4,i2)) + &
                      zb(i2,i1)*zb(i3,i2)*zb(i4,i1)* &
                       (4*mt2 + za(i3,i4)*zb(i4,i3))) - &
                   za(i3,i4)* &
                    (6*za(i3,i4)*zb(i3,i1)**2*zb(i4,i2)**2* &
                       zb(i4,i3) - &
                      zb(i3,i1)*zb(i4,i2)* &
                       (2*zb(i3,i2)*zb(i4,i1)* &
                          (4*mt2 - za(i3,i4)*zb(i4,i3)) + &
                         zb(i2,i1)*zb(i4,i3)* &
                          (8*mt2 + 7*za(i3,i4)*zb(i4,i3))) + &
                      zb(i2,i1)*zb(i4,i3)* &
                       (za(i3,i4)*zb(i2,i1)*zb(i4,i3)**2 + &
                         zb(i3,i2)*zb(i4,i1)* &
                          (4*mt2 + 9*za(i3,i4)*zb(i4,i3))))))) &
          )/(2.*za(i1,i2)**2*za(i3,i4)*zb(i2,i1)**2*zb(i3,i1)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**3* &
          zb(i4,i3))

      return
      end

      function ggHZ_mp_tri1(i1,i2,i3,i4,za,zb,s,mt2)
      implicit none
      complex(8)::ggHZ_mp_tri1,za(1:4,1:4),zb(1:4,1:4)
!---- coefficiient of 2mh box
      integer:: i1,i2,i3,i4
      real(8):: mt2,t,s(1:4,1:4)
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_tri1=(-(za(i1,i2)**3*za(i3,i4)*zb(i2,i1)**2*zb(i4,i1)* &
             zb(i4,i2)) + &
          za(i1,i2)**2*zb(i2,i1)*zb(i4,i1)* &
           (2*za(i1,i3)**2*zb(i2,i1)*zb(i3,i1) + &
             za(i1,i3)*(2*za(i1,i4)*zb(i2,i1)*zb(i4,i1) + &
                (3*za(i2,i4)*zb(i2,i1) - &
                   za(i3,i4)*zb(i3,i1))*zb(i4,i2)) - &
             zb(i4,i2)*(3*za(i1,i4)*za(i2,i3)*zb(i2,i1) + &
                2*za(i3,i4)* &
                 (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))) &
             ) - za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)* &
           zb(i4,i2)*(za(i2,i3)*zb(i3,i2) + &
             za(i2,i4)*zb(i4,i2))*zb(i4,i3) + &
          za(i1,i3)**2*zb(i3,i1)* &
           (za(i2,i4)**2*zb(i4,i1)*zb(i4,i2)**2 - &
             2*za(i2,i3)**2*zb(i2,i1)*zb(i3,i2)*zb(i4,i3) + &
             za(i2,i3)*za(i2,i4)*zb(i4,i2)* &
              (zb(i3,i2)*zb(i4,i1) - 3*zb(i2,i1)*zb(i4,i3))) &
           + za(i1,i3)*(za(i2,i4)*za(i3,i4)*zb(i4,i1)* &
              zb(i4,i2)*(za(i2,i3)*zb(i3,i2) + &
                za(i2,i4)*zb(i4,i2))*zb(i4,i3) - &
             za(i1,i4)*za(i2,i3)* &
              (za(i2,i4)*zb(i4,i1)*zb(i4,i2)* &
                 (zb(i3,i2)*zb(i4,i1) + &
                   3*zb(i2,i1)*zb(i4,i3)) + &
                za(i2,i3)* &
                 (zb(i3,i2)**2*zb(i4,i1)**2 + &
                   2*zb(i2,i1)*zb(i3,i2)*zb(i4,i1)* &
                    zb(i4,i3) - zb(i2,i1)**2*zb(i4,i3)**2))) &
           - za(i1,i2)*(za(i1,i3)**2*zb(i2,i1)*zb(i3,i1)* &
              (-3*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + &
                za(i2,i3)* &
                 (-2*zb(i3,i2)*zb(i4,i1) + &
                   2*zb(i2,i1)*zb(i4,i3))) + &
             zb(i4,i2)*(za(i3,i4)* &
                 (2*za(i2,i3)*zb(i2,i1) + &
                   za(i3,i4)*zb(i4,i1))* &
                 (za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2))* &
                 zb(i4,i3) + &
                2*za(i1,i4)*za(i2,i3)*zb(i2,i1)* &
                 (2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + &
                   za(i2,i3)* &
                    (2*zb(i3,i2)*zb(i4,i1) - &
                      zb(i2,i1)*zb(i4,i3)))) + &
             za(i1,i3)*(za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                 (-2*za(i2,i4)*zb(i4,i1)*zb(i4,i2) + &
                   za(i2,i3)* &
                    (-(zb(i3,i2)*zb(i4,i1)) + &
                      3*zb(i2,i1)*zb(i4,i3))) + &
                zb(i4,i2)* &
                 (za(i2,i4)* &
                    (-4*za(i2,i4)*zb(i2,i1) + &
                      za(i3,i4)*zb(i3,i1))*zb(i4,i1)* &
                    zb(i4,i2) + &
                   za(i2,i3)* &
                    (za(i3,i4)*zb(i3,i1)* &
                       (zb(i3,i2)*zb(i4,i1) - &
                         zb(i2,i1)*zb(i4,i3)) + &
                      2*za(i2,i4)*zb(i2,i1)* &
                       (-2*zb(i3,i2)*zb(i4,i1) + &
                         zb(i2,i1)*zb(i4,i3)))))))/ &
        (2.*s(i3,i4)*za(i1,i2)*zb(i2,i1)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2)
      return
      end

      function ggHZ_mp_tri2(i1,i2,i3,i4,za,zb,s,mt2)
      implicit none
      complex(8):: ggHZ_mp_tri2,za(1:4,1:4),zb(1:4,1:4)
!---- coefficiient of 2mh box
      integer:: i1,i2,i3,i4
      real(8):: mt2,t,s(1:4,1:4)
      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

      ggHZ_mp_tri2= (zb(i4,i1)*(2*za(i1,i3)**3*zb(i2,i1)*zb(i3,i1)* &
             (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) + &
            2*za(i1,i4)*za(i2,i3)*za(i3,i4)*zb(i4,i1)* &
             zb(i4,i2)*(za(i1,i2)*zb(i2,i1) + &
               za(i2,i3)*zb(i3,i2) + za(i2,i4)*zb(i4,i2)) + &
            2*za(i1,i3)*za(i3,i4)* &
             (za(i1,i2)*zb(i2,i1)* &
                (-(za(i1,i4)*zb(i4,i1)**2) + &
                  za(i2,i3)*zb(i3,i1)*zb(i4,i2)) + &
               za(i2,i3)* &
                (zb(i3,i1)*zb(i4,i2)* &
                   (za(i2,i3)*zb(i3,i2) + &
                     za(i2,i4)*zb(i4,i2)) + &
                  za(i1,i4)*zb(i4,i1)* &
                   (-(zb(i3,i2)*zb(i4,i1)) + &
                     zb(i3,i1)*zb(i4,i2)))) + &
            za(i1,i3)**2* &
             (2*za(i1,i4)*zb(i2,i1)*zb(i4,i1)* &
                (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1)) &
                + za(i3,i4)* &
                (-2*za(i1,i2)*zb(i2,i1)*zb(i3,i1)* &
                   zb(i4,i1) + &
                  za(i2,i4)*zb(i4,i1)* &
                   (zb(i3,i2)*zb(i4,i1) - &
                     zb(i3,i1)*zb(i4,i2) + &
                     zb(i2,i1)*zb(i4,i3)) + &
                  za(i2,i3)*zb(i3,i1)* &
                   (-(zb(i3,i2)*zb(i4,i1)) + &
                     zb(i3,i1)*zb(i4,i2) + &
                     zb(i2,i1)*zb(i4,i3))))))/ &
        (2.*za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))**2* &
          zb(i4,i3))

      return
      end


       function ggHZ_mp_3mtri(i1,i2,i3,i4,za,zb,s)
      implicit none
      complex(8):: ggHZ_mp_3mtri,za(1:4,1:4),zb(1:4,1:4)
      integer:: i1,i2,i3,i4
      real(8):: mt2,t,gam
      real(8):: P12(4),P34(4),s(1:4,1:4)
      integer:: nu


      t(i1,i2,i3)=s(i1,i2)+s(i1,i3)+s(i2,i3)

!==== part 1
      ggHZ_mp_3mtri=  (-2*((za(i1,i3)*zb(i3,i1))/2. + &
             (za(i2,i3)*zb(i3,i2))/2. + &
             (za(i1,i4)*zb(i4,i1))/2. + &
             (za(i2,i4)*zb(i4,i2))/2. + &
             Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                  (za(i2,i3)*zb(i3,i2))/2. + &
                  (za(i1,i4)*zb(i4,i1))/2. + &
                  (za(i2,i4)*zb(i4,i2))/2.)**2 - &
               za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))) &
            **2*(1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                zb(i4,i3))/ &
              ((za(i1,i3)*zb(i3,i1))/2. + &
                 (za(i2,i3)*zb(i3,i2))/2. + &
                 (za(i1,i4)*zb(i4,i1))/2. + &
                 (za(i2,i4)*zb(i4,i2))/2. + &
                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                      (za(i2,i3)*zb(i3,i2))/2. + &
                      (za(i1,i4)*zb(i4,i1))/2. + &
                      (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                   za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                    zb(i4,i3)))**2)**8* &
          ((((za(i1,i3)*zb(i3,i1))/2. + &
                 (za(i2,i3)*zb(i3,i2))/2. + &
                 (za(i1,i4)*zb(i4,i1))/2. + &
                 (za(i2,i4)*zb(i4,i2))/2. + &
                 Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                      (za(i2,i3)*zb(i3,i2))/2. + &
                      (za(i1,i4)*zb(i4,i1))/2. + &
                      (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                   za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                    zb(i4,i3)))* &
               (za(i1,i2)*zb(i2,i1) - &
                 (za(i1,i2)*zb(i2,i1)* &
                    (za(i1,i3)*zb(i3,i1) + &
                      za(i1,i4)*zb(i4,i1)))/ &
                  ((za(i1,i3)*zb(i3,i1))/2. + &
                    (za(i2,i3)*zb(i3,i2))/2. + &
                    (za(i1,i4)*zb(i4,i1))/2. + &
                    (za(i2,i4)*zb(i4,i2))/2. + &
                    Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                      za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                       zb(i4,i3))))* &
               (za(i1,i2)*zb(i2,i1)* &
                  (((za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                       (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2)))/ &
                     (1 - &
                        (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2 - &
                    ((za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                       (za(i2,i3)*zb(i3,i2) + &
                         za(i2,i4)*zb(i4,i2) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (1 - &
                        (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2)* &
                  ((za(i1,i2)**2*zb(i2,i1)**2* &
                       (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))**2* &
                       (-(za(i1,i2)*zb(i4,i2)) + &
                         (za(i1,i2)*za(i1,i3)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                       (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2* &
                       (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2* &
                       (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**6) + &
                    ((za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))**2* &
                       (-(za(i2,i3)*zb(i2,i1)) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i1))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                       (za(i1,i2)*zb(i2,i1) - &
                         (za(i1,i2)*zb(i2,i1)* &
                         (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2* &
                       (-(za(i1,i3)*zb(i4,i3)) + &
                         (za(i1,i2)*za(i3,i4)*zb(i4,i2)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (1 - &
                        (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**6) - &
                 (za(i1,i2)*zb(i2,i1)* &
                    (za(i2,i3)*zb(i3,i1) + &
                       za(i2,i4)*zb(i4,i1))**2* &
                    (za(i1,i3)*zb(i3,i1) + &
                      za(i1,i4)*zb(i4,i1) - &
                      (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                       ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                    (za(i1,i2)*zb(i4,i2)* &
                       (((za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                         (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                         (-(za(i2,i3)*zb(i2,i1)) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i1))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i1,i2)*zb(i2,i1) - &
                         (za(i1,i2)*zb(i2,i1)* &
                         (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**4 + &
                         (za(i1,i2)**2*zb(i2,i1)**2* &
                         (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                         (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                         (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2* &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**4)) + &
                      za(i1,i3)*zb(i2,i1)* &
                       (-((za(i1,i2)*zb(i2,i1)* &
                         (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                         (-(za(i1,i2)*zb(i4,i2)) + &
                         (za(i1,i2)*za(i1,i3)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2)/ &
                         (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**4)) + &
                         ((za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                         (za(i1,i2)*zb(i2,i1) - &
                         (za(i1,i2)*zb(i2,i1)* &
                         (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2* &
                         (-(za(i1,i3)*zb(i4,i3)) + &
                         (za(i1,i2)*za(i3,i4)*zb(i4,i2)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**4)))/ &
                  (1 - &
                     (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                        zb(i4,i3))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**3))/ &
             (1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                  zb(i4,i3))/ &
                ((za(i1,i3)*zb(i3,i1))/2. + &
                   (za(i2,i3)*zb(i3,i2))/2. + &
                   (za(i1,i4)*zb(i4,i1))/2. + &
                   (za(i2,i4)*zb(i4,i2))/2. + &
                   Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                        (za(i2,i3)*zb(i3,i2))/2. + &
                        (za(i1,i4)*zb(i4,i1))/2. + &
                        (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                     za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                      zb(i4,i3)))**2) + &
            za(i3,i4)*zb(i4,i3)* &
             (-((za(i1,i2)**4*zb(i2,i1)**4* &
                    (za(i2,i3)*zb(i3,i1) + &
                       za(i2,i4)*zb(i4,i1))**3* &
                    (za(i1,i3)*zb(i3,i2) + &
                      za(i1,i4)*zb(i4,i2))* &
                    (za(i1,i2)*zb(i2,i1) + &
                      (za(i1,i3)*zb(i3,i1))/2. + &
                      (za(i2,i3)*zb(i3,i2))/2. + &
                      (za(i1,i4)*zb(i4,i1))/2. + &
                      (za(i2,i4)*zb(i4,i2))/2. + &
                      Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                        za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                    (-(za(i1,i2)*zb(i4,i2)) + &
                      (za(i1,i2)*za(i1,i3)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                       ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                    (za(i1,i3)*zb(i3,i1) + &
                       za(i1,i4)*zb(i4,i1) - &
                       (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                        ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**3* &
                    (za(i3,i4)*zb(i4,i1) + &
                      (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                       ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                  (((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**4* &
                    (1 - &
                       (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                        ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**9)) + &
               (za(i1,i2)*zb(i2,i1)* &
                  (za(i2,i3)*zb(i3,i1) + &
                     za(i2,i4)*zb(i4,i1))**3* &
                  (-(za(i2,i3)*zb(i2,i1)) - &
                    (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                       zb(i4,i1))/ &
                     ((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                  (za(i1,i2)*zb(i2,i1) - &
                     (za(i1,i2)*zb(i2,i1)* &
                        (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                        (za(i2,i3)*zb(i3,i2))/2. + &
                        (za(i1,i4)*zb(i4,i1))/2. + &
                        (za(i2,i4)*zb(i4,i2))/2. + &
                        Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                        2. + (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2* &
                  (za(i1,i3)*zb(i3,i1) + &
                    za(i1,i4)*zb(i4,i1) - &
                    (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                       zb(i4,i3))/ &
                     ((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                  ((za(i1,i2)**2*zb(i2,i1)*zb(i4,i2)* &
                       (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2)))/ &
                     (1 - &
                       (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                        ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2) - &
                    (za(i1,i2)*zb(i2,i1)* &
                       (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                       (za(i1,i2)*zb(i2,i1) + &
                         (za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                       (-(za(i1,i3)*zb(i4,i3)) + &
                         (za(i1,i2)*za(i3,i4)*zb(i4,i2)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                       (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2)))/ &
                (((za(i1,i3)*zb(i3,i1))/2. + &
                    (za(i2,i3)*zb(i3,i2))/2. + &
                    (za(i1,i4)*zb(i4,i1))/2. + &
                    (za(i2,i4)*zb(i4,i2))/2. + &
                    Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                      za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                       zb(i4,i3)))* &
                  (1 - &
                     (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                        zb(i4,i3))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**7) + &
               ((za(i2,i3)*zb(i3,i1) + &
                     za(i2,i4)*zb(i4,i1))**2* &
                  (za(i1,i2)*zb(i2,i1) - &
                     (za(i1,i2)*zb(i2,i1)* &
                        (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                        (za(i2,i3)*zb(i3,i2))/2. + &
                        (za(i1,i4)*zb(i4,i1))/2. + &
                        (za(i2,i4)*zb(i4,i2))/2. + &
                        Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                        2. + (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**3* &
                  (-(za(i1,i3)*zb(i4,i3)) + &
                    (za(i1,i2)*za(i3,i4)*zb(i4,i2)* &
                       zb(i4,i3))/ &
                     ((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                  ((za(i1,i2)*zb(i2,i1)* &
                       (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                       (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                       (-(za(i2,i3)*zb(i2,i1)) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i1))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (1 - &
                        (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**3 + &
                    ((za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                       (((za(i1,i2)*zb(i2,i1) + &
                         (za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                         (-(za(i2,i3)*zb(i2,i1)) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i1))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i1,i2)*zb(i2,i1) - &
                         (za(i1,i2)*zb(i2,i1)* &
                         (za(i2,i3)*zb(i3,i2) + &
                         za(i2,i4)*zb(i4,i2)))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2 - &
                         za(i1,i2)*zb(i2,i1)* &
                         (-((za(i1,i2)*za(i1,i3)* &
                         zb(i2,i1)**2* &
                         (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1)))/ &
                         (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2))) + &
                         ((-(za(i2,i3)*zb(i2,i1)) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i1))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i2,i3)*zb(i3,i2) + &
                         za(i2,i4)*zb(i4,i2) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2)))/ &
                     (1 - &
                       (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                        ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)))/ &
                (1 - (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                      zb(i4,i3))/ &
                    ((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**6 + &
               (za(i1,i2)**2*zb(i2,i1)**2* &
                  (za(i2,i3)*zb(i3,i1) + &
                     za(i2,i4)*zb(i4,i1))**2* &
                  (za(i1,i2)*zb(i2,i1) - &
                    (za(i1,i2)*zb(i2,i1)* &
                       (za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1)))/ &
                     ((za(i1,i3)*zb(i3,i1))/2. + &
                       (za(i2,i3)*zb(i3,i2))/2. + &
                       (za(i1,i4)*zb(i4,i1))/2. + &
                       (za(i2,i4)*zb(i4,i2))/2. + &
                       Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                  (za(i1,i3)*zb(i3,i1) + &
                     za(i1,i4)*zb(i4,i1) - &
                     (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                        zb(i4,i3))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                        (za(i2,i3)*zb(i3,i2))/2. + &
                        (za(i1,i4)*zb(i4,i1))/2. + &
                        (za(i2,i4)*zb(i4,i2))/2. + &
                        Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                        2. + (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))**2* &
                  ((za(i1,i2)**3*zb(i2,i1)**2* &
                       (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))*zb(i4,i2)* &
                       (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                       (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                     (((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                       (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**3) + &
                    ((-(za(i1,i2)*zb(i4,i2)) + &
                         (za(i1,i2)*za(i1,i3)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                       ((za(i1,i2)*zb(i2,i1)* &
                         (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1))* &
                         (za(i1,i3)*zb(i3,i2) + &
                         za(i1,i4)*zb(i4,i2))* &
                         (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**3 + &
                         ((za(i1,i3)*zb(i3,i1) + &
                         za(i1,i4)*zb(i4,i1) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (((za(i1,i2)*zb(i2,i1) + &
                         (za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))* &
                         (za(i1,i2)*zb(i2,i1) - &
                         (za(i1,i2)*zb(i2,i1)* &
                         (za(i2,i3)*zb(i3,i2) + &
                         za(i2,i4)*zb(i4,i2)))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2 - &
                         za(i1,i2)*zb(i2,i1)* &
                         ((za(i1,i3)*zb(i2,i1)* &
                         (za(i2,i3)*zb(i3,i1) + &
                         za(i2,i4)*zb(i4,i1)))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2) + &
                         ((za(i2,i3)*zb(i3,i2) + &
                         za(i2,i4)*zb(i4,i2) - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))))* &
                         (za(i3,i4)*zb(i4,i1) + &
                         (za(i2,i3)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**2)))/ &
                         (1 - &
                         (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                         ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)))/ &
                     (1 - &
                       (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3))/ &
                        ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)))/ &
                (((za(i1,i3)*zb(i3,i1))/2. + &
                     (za(i2,i3)*zb(i3,i2))/2. + &
                     (za(i1,i4)*zb(i4,i1))/2. + &
                     (za(i2,i4)*zb(i4,i2))/2. + &
                     Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                       za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                        zb(i4,i3)))**2* &
                  (1 - &
                     (za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                        zb(i4,i3))/ &
                      ((za(i1,i3)*zb(i3,i1))/2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2. + &
                         Sqrt(((za(i1,i3)*zb(i3,i1))/ &
                         2. + &
                         (za(i2,i3)*zb(i3,i2))/2. + &
                         (za(i1,i4)*zb(i4,i1))/2. + &
                         (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                         za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                         zb(i4,i3)))**2)**5))))/ &
        (za(i1,i2)**2*zb(i2,i1)**2* &
          (za(i2,i3)*zb(i3,i1) + za(i2,i4)*zb(i4,i1))** &
           4*(za(i1,i2)*zb(i2,i1) - &
             (za(i1,i2)*zb(i2,i1)* &
                (za(i1,i3)*zb(i3,i1) + &
                  za(i1,i4)*zb(i4,i1)))/ &
              ((za(i1,i3)*zb(i3,i1))/2. + &
                (za(i2,i3)*zb(i3,i2))/2. + &
                (za(i1,i4)*zb(i4,i1))/2. + &
                (za(i2,i4)*zb(i4,i2))/2. + &
                Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                     (za(i2,i3)*zb(i3,i2))/2. + &
                     (za(i1,i4)*zb(i4,i1))/2. + &
                     (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                  za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                   zb(i4,i3))))**2* &
          (za(i1,i3)*zb(i3,i1) + za(i1,i4)*zb(i4,i1) - &
             (za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3))/ &
              ((za(i1,i3)*zb(i3,i1))/2. + &
                (za(i2,i3)*zb(i3,i2))/2. + &
                (za(i1,i4)*zb(i4,i1))/2. + &
                (za(i2,i4)*zb(i4,i2))/2. + &
                Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                     (za(i2,i3)*zb(i3,i2))/2. + &
                     (za(i1,i4)*zb(i4,i1))/2. + &
                     (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                  za(i1,i2)*za(i3,i4)*zb(i2,i1)* &
                   zb(i4,i3))))**2* &
          (-(za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3)) + &
            ((za(i1,i3)*zb(i3,i1))/2. + &
               (za(i2,i3)*zb(i3,i2))/2. + &
               (za(i1,i4)*zb(i4,i1))/2. + &
               (za(i2,i4)*zb(i4,i2))/2. + &
               Sqrt(((za(i1,i3)*zb(i3,i1))/2. + &
                    (za(i2,i3)*zb(i3,i2))/2. + &
                    (za(i1,i4)*zb(i4,i1))/2. + &
                    (za(i2,i4)*zb(i4,i2))/2.)**2 - &
                 za(i1,i2)*za(i3,i4)*zb(i2,i1)*zb(i4,i3) &
                 ))**2))

!+==== part 2
!=====normalization
      ggHZ_mp_3mtri=-ggHZ_mp_3mtri/(two*s(i3,i4))

      return
      end
! ------------------ END MCFM SM check  ------------      






end module ModVHiggs
!!--YaofuZhou-----------------------------------------
