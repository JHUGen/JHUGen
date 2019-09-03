!--YaofuZhou-----------------------------------------
module ModVHiggs
  use ModParameters
  use ModMisc
  implicit none
  private


  !----- notation for subroutines
  public :: EvalAmp_VHiggs

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

         if (includeVprime) then
            ghzzp1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q3_q3,q5_q5)
            ghzzp2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q3_q3,q5_q5)
            ghzzp3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q3_q3,q5_q5)
            ghzzp4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q3_q3,q5_q5)
         endif
      else !if(useA(2)) then
         ghz1_dyn = HVVSpinZeroDynamicCoupling(5,0d0,q4_q4,q5_q5)
         ghz2_dyn = HVVSpinZeroDynamicCoupling(6,0d0,q4_q4,q5_q5)
         ghz3_dyn = HVVSpinZeroDynamicCoupling(7,0d0,q4_q4,q5_q5)
         ghz4_dyn = HVVSpinZeroDynamicCoupling(8,0d0,q4_q4,q5_q5)

         if (includeVprime) then
            ghzpz1_dyn = HVVSpinZeroDynamicCoupling(20,0d0,q4_q4,q5_q5)
            ghzpz2_dyn = HVVSpinZeroDynamicCoupling(21,0d0,q4_q4,q5_q5)
            ghzpz3_dyn = HVVSpinZeroDynamicCoupling(22,0d0,q4_q4,q5_q5)
            ghzpz4_dyn = HVVSpinZeroDynamicCoupling(23,0d0,q4_q4,q5_q5)
         endif
      endif

      gVVS1 = ghz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda**2 )
      gVVS2 = -( 2d0*ghz2_dyn + ghz3_dyn*qq/Lambda**2 )
      gVVP = -2d0*ghz4_dyn

      if (includeVprime) then
         if(.not.useA(1) .and. .not.useA(2)) then
            gVVpS1 = ghzzp1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpS2 = -( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpP = -2d0*ghzzp4_dyn

            gVpVS1 = ghzpz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVS2 = -( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVP = -2d0*ghzpz4_dyn

            gVpVpS1 = ghzpzp1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzpzp2_dyn + ghzpzp3_dyn*qq/Lambda**2 )
            gVpVpS2 = -( 2d0*ghzpzp2_dyn + ghzpzp3_dyn*qq/Lambda**2 )
            gVpVpP = -2d0*ghzpzp4_dyn
         else if (useA(1)) then
            gVVpS1 = ghzzp1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpS2 = -( 2d0*ghzzp2_dyn + ghzzp3_dyn*qq/Lambda**2 )
            gVVpP = -2d0*ghzzp4_dyn
         else! if (useA(2)) then
            gVpVS1 = ghzpz1_dyn*(mass(3,1)**2) + qq * ( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVS2 = -( 2d0*ghzpz2_dyn + ghzpz3_dyn*qq/Lambda**2 )
            gVpVP = -2d0*ghzpz4_dyn
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




end module ModVHiggs
!!--YaofuZhou-----------------------------------------
