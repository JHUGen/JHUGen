      !-- H -> [l^-(i3) lb(i4)] [l^-(i5) lb(i6)]
      function anomhzzamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'mxpart.f'
      include 'masses.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhzzamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4,LambdaJHBSM
      double complex ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3

c------ HZZ DECAY CONVENTIONS
      IF( AllowAnomalousCouplings.eq.1 ) THEN

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(i3,i5)*zb(i4,i6)*shat
      FFa2 = -0.5d0*za(i3,i5)**2*zb(i3,i6)*zb(i4,i5)
     &       -0.5d0*za(i3,i5)*za(i3,i6)*zb(i3,i6)*zb(i4,i6)
     &       -0.5d0*za(i3,i5)*za(i4,i5)*zb(i4,i5)*zb(i4,i6)
     &       -0.5d0*za(i3,i6)*za(i4,i5)*zb(i4,i6)**2
      FFa3 =  0.5d0*za(i3,i4)*za(i5,i6)*zb(i4,i6)**2
     &       -0.5d0*za(i3,i5)**2*zb(i3,i4)*zb(i5,i6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

c--- q^2 dependent couplings
      call HVVSpinZeroDynCoupl(ghz1_dyn,1,jh,q3_q3,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz2_dyn,2,jh,q3_q3,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz3_dyn,3,jh,q3_q3,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz4_dyn,4,jh,q3_q3,q4_q4,shat,.false.)

      if(jh.eq.1) then
      LambdaJHBSM=LambdaBSM
      else
      LambdaJHBSM=Lambda2BSM
      endif

      aa1 =ghz1_dyn*zmass**2/shat
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaJHBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaJHBSM**2
      aa3 =-2d0*ghz4_dyn


      aa1 = aa1 / zmass**2 !-- F
      aa2 = aa2 / zmass**2 !-- F
      aa3 = aa3 / zmass**2 !-- F

      anomhzzamp = ( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )

      ELSE

      anomhzzamp = za(i3,i5)*zb(i4,i6)

      ENDIF

      return

      end


      !-- H -> [l^-(i3) lb(i4)] [(gamma* ->) l^-(i5) lb(i6)]
      function anomhzaamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhzaamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4,LambdaJHBSM
      double complex ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3
      logical checkHZAcouplings

      anomhzaamp = czip ! No ZA->4f amplitude at LO SM

c------ HZZ DECAY CONVENTIONS
      IF(AllowAnomalousCouplings.eq.1) THEN
      if(checkHZAcouplings(jh)) then

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(i3,i5)*zb(i4,i6)*shat
      FFa2 = -0.5d0*za(i3,i5)**2*zb(i3,i6)*zb(i4,i5)
     &       -0.5d0*za(i3,i5)*za(i3,i6)*zb(i3,i6)*zb(i4,i6)
     &       -0.5d0*za(i3,i5)*za(i4,i5)*zb(i4,i5)*zb(i4,i6)
     &       -0.5d0*za(i3,i6)*za(i4,i5)*zb(i4,i6)**2
      FFa3 =  0.5d0*za(i3,i4)*za(i5,i6)*zb(i4,i6)**2
     &       -0.5d0*za(i3,i5)**2*zb(i3,i4)*zb(i5,i6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

c--- q^2 dependent couplings
      ! q4_q4==q**2 of gamma*
      call HVVSpinZeroDynCoupl(ghz1_dyn,5,jh,0d0,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz2_dyn,6,jh,0d0,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz3_dyn,7,jh,0d0,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz4_dyn,8,jh,0d0,q4_q4,shat,.false.)

      if(jh.eq.1) then
      LambdaJHBSM=LambdaBSM
      else
      LambdaJHBSM=Lambda2BSM
      endif

      aa1 =ghz1_dyn*zmass**2/shat
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaJHBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaJHBSM**2
      aa3 =-2d0*ghz4_dyn


      aa1 = aa1 / zmass**2 !-- F
      aa2 = aa2 / zmass**2 !-- F
      aa3 = aa3 / zmass**2 !-- F

      anomhzaamp = ( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )
      endif
      ENDIF

      return

      end


      !-- H -> [(gamma* ->) l^-(i3) lb(i4)] [(gamma* ->) l^-(i5) lb(i6)]
      function anomhaaamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhaaamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4,LambdaJHBSM
      double complex ghz2_dyn,ghz3_dyn,ghz4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3
      logical checkHAAcouplings

      anomhaaamp = czip ! No AA->4f amplitude at LO SM

c------ HZZ DECAY CONVENTIONS
      IF(AllowAnomalousCouplings.eq.1) THEN
      if(checkHAAcouplings(jh)) then

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(i3,i5)*zb(i4,i6)*shat
      FFa2 = -0.5d0*za(i3,i5)**2*zb(i3,i6)*zb(i4,i5)
     &       -0.5d0*za(i3,i5)*za(i3,i6)*zb(i3,i6)*zb(i4,i6)
     &       -0.5d0*za(i3,i5)*za(i4,i5)*zb(i4,i5)*zb(i4,i6)
     &       -0.5d0*za(i3,i6)*za(i4,i5)*zb(i4,i6)**2
      FFa3 =  0.5d0*za(i3,i4)*za(i5,i6)*zb(i4,i6)**2
     &       -0.5d0*za(i3,i5)**2*zb(i3,i4)*zb(i5,i6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

c--- q^2 dependence in couplings does not exist in gsgs
      call HVVSpinZeroDynCoupl(ghz2_dyn,9,jh,q3_q3,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz3_dyn,10,jh,q3_q3,q4_q4,shat,.false.)
      call HVVSpinZeroDynCoupl(ghz4_dyn,11,jh,q3_q3,q4_q4,shat,.false.)

      if(jh.eq.1) then
      LambdaJHBSM=LambdaBSM
      else
      LambdaJHBSM=Lambda2BSM
      endif

      aa1 = !ghz1_dyn*zmass**2/shat ! No ghz1
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaJHBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaJHBSM**2
      aa3 =-2d0*ghz4_dyn


      ! Divide by zmass**2 just like in ZZ or ZA
      ! although it has no meaning here
      aa1 = aa1 / zmass**2 !-- F
      aa2 = aa2 / zmass**2 !-- F
      aa3 = aa3 / zmass**2 !-- F

      anomhaaamp = ( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )

      endif
      ENDIF

      return

      end


      !-- H -> [nu(i3) l^+(i4)] [l^-(i5) nub(i6)]
      function anomhwwamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhwwamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4,LambdaJHBSM
      double complex ghw1_dyn,ghw2_dyn,ghw3_dyn,ghw4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3

c------ HZZ=-HWW DECAY CONVENTIONS
      IF( AllowAnomalousCouplings.eq.1 ) THEN

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(i3,i5)*zb(i4,i6)*shat
      FFa2 = -0.5d0*za(i3,i5)**2*zb(i3,i6)*zb(i4,i5)
     &       -0.5d0*za(i3,i5)*za(i3,i6)*zb(i3,i6)*zb(i4,i6)
     &       -0.5d0*za(i3,i5)*za(i4,i5)*zb(i4,i5)*zb(i4,i6)
     &       -0.5d0*za(i3,i6)*za(i4,i5)*zb(i4,i6)**2
      FFa3 =  0.5d0*za(i3,i4)*za(i5,i6)*zb(i4,i6)**2
     &       -0.5d0*za(i3,i5)**2*zb(i3,i4)*zb(i5,i6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

c--- q^2-dependent couplings
      call HVVSpinZeroDynCoupl(ghw1_dyn,1,jh,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw2_dyn,2,jh,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw3_dyn,3,jh,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw4_dyn,4,jh,q3_q3,q4_q4,shat,.true.)

      if(jh.eq.1) then
      LambdaJHBSM=LambdaBSM
      else
      LambdaJHBSM=Lambda2BSM
      endif

      aa1 =ghw1_dyn*wmass**2/shat
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghw2_dyn
     &       +ghw3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaJHBSM**2)
      aa2 =-2d0*ghw2_dyn
     &     -ghw3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaJHBSM**2
      aa3 =-2d0*ghw4_dyn

      aa1 = aa1 / wmass**2
      aa2 = aa2 / wmass**2
      aa3 = aa3 / wmass**2

      anomhwwamp = ( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )
      !print *,"anomhwwamp(",i3,i4,i5,i6,")=",anomhwwamp

      ELSE

      anomhwwamp = za(i3,i5)*zb(i4,i6)

      ENDIF

      return

      end


      !-- g(i1) g(i2) =|>- + >- H
      subroutine anomhggvtxamp(i1,i2,jh,za,zb,hggvtxamp)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer i1,i2,jh ! jh is Higgs 1, 2
      integer iq,igen,iv
      double precision mfsq(2,2) ! (qgen,b/t)
      double precision tauinv(2,2) ! (qgen,b/t)
      double complex a1(2),a3(2) ! (qgen)
      double complex kappaj(2,2) ! (qgen,b/t)
      double complex kappaj_tilde(2,2) ! (qgen,b/t)
      double complex qlI3,C0mXq
      double complex hggvtxamp(2,2,2),hggpointvtx(2,2)

      mfsq(1,1)=mb**2
      mfsq(1,2)=mt**2
      mfsq(2,1)=mb_4gen**2
      mfsq(2,2)=mt_4gen**2
      tauinv(:,:)=4d0*mfsq(:,:)/s(i1,i2)

      hggvtxamp(:,:,:)=czip
      hggpointvtx(:,:)=czip
      a1(:)=czip
      a3(:)=czip
      kappaj(:,:)=czip
      kappaj_tilde(:,:)=czip

      IF( AllowAnomalousCouplings.eq.1 ) THEN

      if(jh.eq.1) then

      a1(1) = ghg2+ghg3*s(i1,i2)/4d0/LambdaBSM**2
      a3(1) = -2d0*ghg4
      a1(2)= ghg2_4gen+ghg3_4gen*s(i1,i2)/4d0/LambdaBSM**2
      a3(2)= -2d0*ghg4_4gen

      kappaj(1,1) = kappa_bot
      kappaj_tilde(1,1) = kappa_tilde_bot
      kappaj(1,2) = kappa_top
      kappaj_tilde(1,2) = kappa_tilde_top
      kappaj(2,1) = kappa_4gen_bot
      kappaj_tilde(2,1) = kappa_tilde_4gen_bot
      kappaj(2,2) = kappa_4gen_top
      kappaj_tilde(2,2) = kappa_tilde_4gen_top

      else

      a1(1) = gh2g2+gh2g3*s(i1,i2)/4d0/Lambda2BSM**2
      a3(1) = -2d0*gh2g4
      a1(2)= gh2g2_4gen+gh2g3_4gen*s(i1,i2)/4d0/Lambda2BSM**2
      a3(2)= -2d0*gh2g4_4gen

      kappaj(1,1) = kappa2_bot
      kappaj_tilde(1,1) = kappa2_tilde_bot
      kappaj(1,2) = kappa2_top
      kappaj_tilde(1,2) = kappa2_tilde_top
      kappaj(2,1) = kappa2_4gen_bot
      kappaj_tilde(2,1) = kappa2_tilde_4gen_bot
      kappaj(2,2) = kappa2_4gen_top
      kappaj_tilde(2,2) = kappa2_tilde_4gen_top

      endif

      ELSE

      kappaj(1,1) = cone
      kappaj(1,2) = cone

      ENDIF

      do igen=1,2
      ! Compute the point interaction
      hggpointvtx(1,1) =
     & s(i1,i2)*(
     & a1(igen)
     & - a3(igen)*0.5d0*im*za(i1,i2)*zb(i1,i2)/s(i1,i2)
     & )/3d0
      hggpointvtx(2,2) =
     & s(i1,i2)*(
     & a1(igen)
     & + a3(igen)*0.5d0*im*za(i1,i2)*zb(i1,i2)/s(i1,i2)
     & )/3d0

      do iq=1,2
      iv = mod((iq-1),2)+1

      ! Add one point interaction into the SM-top loop
      if (iq.eq.2) then
      hggvtxamp(iv,:,:) = hggvtxamp(iv,:,:) +
     & hggpointvtx(:,:)
      endif

      C0mXq =
     & qlI3(
     & zip,zip,s(i1,i2),
     & mfsq(igen,iq),mfsq(igen,iq),mfsq(igen,iq),
     & musq,0
     & )

      hggvtxamp(iv,1,1) = hggvtxamp(iv,1,1) +
     & mfsq(igen,iq)*(
     &  kappaj(igen,iq)*(
     &    2d0-C0mXq*s(i1,i2)*(1d0-tauinv(igen,iq))
     &  )
     &  -kappaj_tilde(igen,iq)*(-im)*C0mXq*s(i1,i2)
     & )

      hggvtxamp(iv,2,2) = hggvtxamp(iv,2,2) +
     & mfsq(igen,iq)*(
     &  kappaj(igen,iq)*(
     &    2d0-C0mXq*s(i1,i2)*(1d0-tauinv(igen,iq))
     &  )
     &  +kappaj_tilde(igen,iq)*(-im)*C0mXq*s(i1,i2)
     & )

      enddo
      enddo

      ! Multiply by the tensor structure
      do iv=1,2
      hggvtxamp(iv,1,1)=(hggvtxamp(iv,1,1))*za(i1,i2)/zb(i1,i2) 
      hggvtxamp(iv,2,2)=(hggvtxamp(iv,2,2))*zb(i1,i2)/za(i1,i2)
      enddo

      return
      end

c--- Interpolation function for form factors of trilinear
      subroutine lin_interpolate(xvals, yvals, length, x, y)
      implicit none
      integer :: length, i, t
      real*8 :: xvals(length), yvals(length), x0, x1, y0, y1
      real*8, intent(in) :: x
      real*8, intent(out) :: y
      if (x<xvals(1)) then
         y=0
      end if
      if (x>xvals(size(xvals))) then
         y=0
      end if
         
      do i = 1, size(xvals)
         if (xvals(i) > x) then
            x0 = xvals(i-1)
            x1 = xvals(i)
            y0 = yvals(i-1)
            y1 = yvals(i)
            exit
         end if
      end do
        
      if (x>xvals(1)) then
         if (x<xvals(size(xvals))) then
            y = (y1-y0)/(x1-x0)*(x-x0) + y0
         end if
      end if  

      return 
      end subroutine lin_interpolate

      subroutine anomhggvtxamp_c6(za,zb,ggHmt_c6)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'AnomTriLinear.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcouple.f'
      double complex hggvtxamp(2,2,2)
      double complex ggHmt_c6(2,2)
      double precision mb2,mt2,sinthw
      double precision qhsq, q1sq, q2sq, MH2, MZ2, dB0h, dZh, Z
      double precision Fre_x(314), Fre_y(314), Fim_x(314), Fim_y(314)
      double precision Fre, Fim, F_1l, F_2l, sqrts 
      double complex qlI3,C0mt
c--- squared masses and sin(thetaw)     
      mt2=mt**2
      mb2=mb**2
      MH2  = hmass**2
      sinthw=sqrt(xw)

c--- Needed factor 
      dB0h = (-9 + 2*Sqrt(3.)*Pi)/(9.*MH2)


c--- Amplitudes for production 
      C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)

c--- c6 correction to ggH vertex      

c--- The real part of the 2-loop form factor 
      
      Fre_x  =(/ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0,
     &12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 
     &24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0,
     &38.0, 40.0, 42.0, 44.0, 46.0, 48.0,  
     &50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0,
     &64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 
     &76.0, 78.0, 80.0, 82.0, 84.0, 86.0, 88.0,
     &90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 
     &102.0, 104.0, 106.0, 108.0, 110.0, 
     $112.0, 114.0, 116.0, 118.0, 120.0, 122.0, 
     &124.0, 126.0, 128.0, 130.0, 132.0, 
     $134.0, 136.0, 138.0, 140.0, 142.0, 144.0, 
     &146.0, 148.0, 150.0, 152.0, 154.0, 
     $156.0, 160.0, 170.0, 180.0, 190.0, 200.0,
     &210.0, 220.0, 230.0, 240.0, 251.0, 
     $261.0, 271.0, 281.0, 291.0, 301.0, 311.0, 
     &321.0, 331.0, 341.0, 351.0, 361.0, 
     $371.0, 381.0, 391.0, 401.0, 411.0, 421.0, 
     &431.0, 441.0, 451.0, 461.0, 471.0, 
     $481.0, 491.0, 501.0, 511.0, 521.0, 531.0, 
     &541.0, 551.0, 561.0, 571.0, 581.0, 
     $591.0, 601.0, 611.0, 621.0, 631.0, 641.0, 
     &651.0, 661.0, 671.0, 681.0, 691.0, 
     $701.0, 711.0, 721.0, 731.0, 741.0, 751.0, 
     &761.0, 771.0, 781.0, 791.0, 801.0, 
     $811.0, 821.0, 831.0, 841.0, 851.0, 861.0, 
     &871.0, 881.0, 891.0, 901.0, 911.0, 
     $921.0, 931.0, 941.0, 951.0, 961.0, 971.0, 
     &981.0, 991.0, 1001.0, 1011.0, 1021.0,
     &1031.0, 1041.0, 1051.0, 1061.0, 1071.0, 
     &1081.0, 1091.0, 1101.0, 1111.0, 1121.0, 
     &1131.0, 1141.0, 1151.0, 1161.0, 1171.0, 
     &1181.0, 1191.0, 1201.0, 1211.0, 1221.0, 
     &1231.0, 1241.0, 1251.0, 1261.0, 1271.0, 
     &1281.0, 1291.0, 1301.0, 1311.0, 1321.0, 
     &1331.0, 1341.0, 1351.0, 1361.0, 1371.0, 
     &1381.0, 1391.0, 1401.0, 1411.0, 1421.0, 
     &1431.0, 1441.0, 1451.0, 1461.0, 1471.0, 
     &1481.0, 1491.0, 1501.0, 1511.0, 1521.0, 
     &1531.0, 1541.0, 1551.0, 1561.0, 1571.0, 
     &1581.0, 1591.0, 1601.0, 1611.0, 1621.0, 
     &1631.0, 1641.0, 1651.0, 1661.0, 1671.0, 
     &1681.0, 1691.0, 1701.0, 1711.0, 1721.0, 
     &1731.0, 1741.0, 1751.0, 1761.0, 1771.0, 
     &1781.0, 1791.0, 1801.0, 1811.0, 1821.0, 
     &1831.0, 1841.0, 1851.0, 1861.0, 1871.0, 
     &1881.0, 1891.0, 1901.0, 1911.0, 1921.0, 
     &1931.0, 1941.0, 1951.0, 1961.0, 1971.0, 
     &1981.0, 1991.0, 2001.0, 2011.0, 2021.0, 
     &2031.0, 2041.0, 2051.0, 2061.0, 2071.0, 
     &2081.0, 2091.0, 2101.0, 2111.0, 2121.0, 
     &2131.0, 2141.0, 2151.0, 2161.0, 2171.0, 
     &2181.0, 2191.0, 2201.0, 2211.0, 2221.0, 
     &2231.0, 2241.0, 2251.0, 2261.0, 2271.0, 
     &2281.0, 2291.0, 2301.0, 2311.0, 2321.0, 
     &2331.0, 2341.0, 2351.0, 2361.0, 2371.0, 
     &2381.0, 2391.0, 2401.0, 2411.0, 2421.0, 
     &2431.0, 2441.0, 2451.0, 2461.0, 2471.0, 
     &2481.0, 2491.0, 2501.0 /)

      Fre_y   = 
     &[ 1.2078594486413152, 1.2078496632938285  , 1.2078204693881747 , 
     &1.2077723533243454  , 1.2077061257360673  ,                    
     &1.2076229214478087  , 1.2075241994273913 , 1.2074117427469975 ,
     &1.2072876585690329  , 1.207154378176958   ,                    
     &1.20701465707486    , 1.206871575183256  , 1.2067285371622782 ,
     &1.2065892728971799  , 1.2064578381848559  ,                    
     &1.2063386156639953  , 1.2062363160353893 , 1.2061559796230659 ,
     &1.2061029783311261  , 1.2060830180556623  ,                    
     &1.2061021416157895  , 1.20616673227288   , 1.2062835179124352 ,
     &1.2064595759688503  , 1.2067023391797087  ,                    
     &1.2070196022631496  , 1.2074195296196466 , 1.2079106641680006 ,
     &1.2085019374350474  , 1.2092026810292669  ,                    
     &1.2100226396406242  , 1.2109719857228145 , 1.2120613360295622 ,
     &1.2133017701945656  , 1.2147048515648955  ,                    
     &1.2162826505210171  , 1.218047770543305  , 1.2200133773156998 ,
     &1.22219323119265    , 1.2246017233964892  ,                    
     &1.2272539163599994  , 1.2301655886840799 , 1.2333532852448503 ,
     &1.236834373059605   , 1.2406271036090377  ,                    
     &1.2447506824165067  , 1.2492253468068604 , 1.2540724529112839 ,
     &1.259314573155359   , 1.2649756056707862  ,                    
     &1.2710808973141638  , 1.277657382267782  , 1.2847337385487305 ,
     &1.2923405651781432  , 1.300510583280025   ,                    
     &1.3092788650120935  , 1.3186830950094082 , 1.3287638699838915 ,
     &1.3395650433201616  , 1.3511341230071356  ,                    
     &1.363522733134738   , 1.3767871515849026 , 1.39098893961737 ,
     &1.4061956830140412  , 1.422481869605585   ,                    
     &1.43992993478623    , 1.4586315156291678 , 1.4786889663083291 ,
     &1.500217203966295   , 1.5233459767931212  ,                    
     &1.5482226776745822  , 1.5750158715762972 , 1.6039197694738974 ,
     &1.635159976675711   , 1.6690009860616684  ,                    
     &1.7057561060390616  , 1.745800859020151  , 1.7895914487761762 ,
     &1.8376908415048878  , 1.9479279734978134  ,                    
     &2.2688816979556266  , 2.6392072429694693 , 3.0667545765855952 ,
     &3.5788341367470506  , 4.192979580691294   ,                    
     &4.9519653618119825  , 5.920999220182388  , 7.281994382695284 ,
     &10.759542275736175  , 11.078208363391356  ,                    
     &11.410183974230948  , 11.672515537835944  , 11.914493211950735 , 
     &12.177126915017174  , 12.463166922896807  ,                    
     &12.738478654069517  , 13.031788523125512  , 13.395904295649041 , 
     &13.358921050601348  , 11.180243634688116  ,                    
     &8.013795654177743   , 4.388898769280723   , 0.6002810584234002 , 
     &-3.180361565005783  , -6.843642433699764  ,                    
     &-10.32029963673472  , -13.570501680239873 , -16.57100288221112 , 
     &-19.316115583372405 , -21.800498188125328 ,                    
     &-24.02524693979271  , -25.994853757057605 , -27.714889085270475, 
     &-29.205650551870814 , -30.48121156584683  ,                    
     &-31.554820658890886 , -32.45086799018988  , -33.190498970308276, 
     &-33.785073898623565 , -34.2476927611582   ,                    
     &-34.59379150145761  , -34.833240257729116 , -34.97849340070205 , 
     &-35.04052825439196  , -35.028852281082365 ,                    
     &-34.94978931002431  , -34.811578318055965 , -34.62022227815581 , 
     &-34.38164606463548  , -34.10232702598602  ,                    
     &-33.78639802399073  , -33.43904866006925  , -33.06381367295392 , 
     &-32.66384003585391  , -32.242548185473446 ,                    
     &-31.801815757223203 , -31.345841835557977 , -30.876818439647252, 
     &-30.395552802224007 , -29.906356285456727 ,                    
     &-29.408286475203667 , -28.903962548412018 , -28.394503199531556, 
     &-27.88237936919367  , -27.367569248972245 ,                    
     &-26.851387866447816 , -26.334884717267666 , -25.818588011227916, 
     &-25.303478593899793 , -24.788936379653933 ,                    
     &-24.276874810370753 , -23.767530987363646 , -23.260685679155923, 
     &-22.75732421110552  , -22.25813534794967  ,                    
     &-21.763468614121056 , -21.272219240452625 , -20.785908226066674, 
     &-20.30378746366105  , -19.826662408423044 ,                    
     &-19.355649387030613 , -18.88925122690857  , -18.42798662838671 , 
     &-17.972507507488764 , -17.522148776860202 ,                    
     &-17.077004713089007 , -16.637768778970422 , -16.204203481408207, 
     &-15.77629661383223  , -15.354504410178608 ,                    
     &-14.93838917246009  , -14.527332589843287 , -14.122076667485034, 
     &-13.723059708842758 , -13.328859805136675 ,                    
     &-12.94101668769741  , -12.557560427415934 , -12.180189712437757, 
     &-11.807675901475747 , -11.441321495014545 ,                    
     &-11.079517461214158 , -10.723179074799194 , -10.371822579388668, 
     &-10.02524400358767  , -9.684091205315918  ,                    
     &-9.348399462200453  , -9.016925131409808  , -8.690916174744588 , 
     &-8.368555262099479  , -8.051298659162383  ,                    
     &-7.739315948307846  , -7.431344320156965  , -7.127521952026591 , 
     &-6.828513203862644  , -6.533799803037556  ,                    
     &-6.243756108336018  , -5.957402098841168  , -5.675127211203429 , 
     &-5.396685136256635  , -5.123768401597527  ,                    
     &-4.85394954478054   , -4.587949259856917  , -4.324203971313953 , 
     &-4.066321298277019  , -3.811492902293355  ,                    
     &-3.5610386317334624 , -3.312320708580716  , -3.069707520950964 , 
     &-2.829086099899587  , -2.591576546496253  ,                    
     &-2.358080844926375  , -2.1274313088386894 , -1.9009174203892654, 
     &-1.6769547828277562 , -1.4562724023045492 ,                    
     &-1.2385090559505734 , -1.0235315200021113 , -0.8124005670154241, 
     &-0.6038016350498732 , -0.3972151697451124 ,                    
     &-0.1945931038893085 , 0.006784394930668096, 0.20408990097333213, 
     &0.3984838280334188  , 0.5902104557188247  ,                    
     &0.7800204956313667  , 0.9667939068706186  , 1.1508788726377015 , 
     &1.3335451347600866  , 1.5116852614705305  ,                    
     &1.6883307551718434  , 1.8614122268390232  , 2.033408027862708  , 
     &2.2033657601055485  , 2.3702829142351587  ,                    
     &2.536663733306287   , 2.698317496636974   , 2.8591790207209082 , 
     &3.0170262402970884  , 3.1740602070150503  ,                    
     &3.32783142770536    , 3.4790363536382864  , 3.630917684096081  , 
     &3.7781757616593694  , 3.9250133785467662  ,                    
     &4.069491181116724   , 4.210025471178852   , 4.35144734965057   , 
     &4.4891050727896555  , 4.626764436352352   ,                    
     &4.760863062370292   , 4.895561297511757   , 5.025242477516087  , 
     &5.156993675234493   , 5.282318187258721   ,                    
     &5.409237780603455   , 5.533557082784054   , 5.654886573786325  , 
     &5.776949192015766   , 5.89861988018190    ,                    
     &6.014801642823237   , 6.132857247798013   , 6.248389653590554  , 
     &6.363490038934817   , 6.477104972780086   ,                    
     &6.58656839925101    , 6.697934388458292   , 6.806827833046982  , 
     &6.912142688984153   , 7.019671790336791   ,                    
     &7.1250715166335885  , 7.230252304320101   , 7.332913649534926  , 
     &7.4343108631583705  , 7.534297868760447   ,                    
     &7.632680499106008   , 7.731782431197287   , 7.828422270802416  , 
     &7.922653350286592   , 8.018183788966036   ,                    
     &8.111729108023441   , 8.202154021260538   , 8.294502553101376  , 
     &8.384128893089471   , 8.472977927501615   ,                    
     &8.563662419392331   , 8.649260306123347   , 8.737915950630759  , 
     &8.820754800244666   , 8.907254199034474   ,                    
     &8.989149046721737   , 9.07107755789215    , 9.15475869312462   , 
     &9.235001084273023   , 9.313362822905678   ,                    
     &9.391042434246772   , 9.471087979878535   , 9.547786305549456  , 
     &9.62527678662758    , 9.699686094446976   ,                    
     &9.774894483690897   , 9.851018773059625   , 9.921578503104827  , 
     &9.991144252128688   , 10.065651267547478  ,                    
     &10.137883647975622  , 10.20703025083462   , 10.278716711904785 , 
     &10.344828135304034  , 10.415921322713487  ,                    
     &10.483065816359064  , 10.548154969712368  , 10.61414302801367  , 
     &10.681484422285573 ]

c--- The imaginary part of the 2-loop form factor       

      Fim_x  = 
     &[ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0,
     &12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 
     &24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0,
     &38.0, 40.0, 42.0, 44.0, 46.0, 48.0,  
     &50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0,
     &64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 
     &76.0, 78.0, 80.0, 82.0, 84.0, 86.0, 88.0,
     &90.0, 92.0, 94.0, 96.0, 98.0, 100.0, 
     &102.0, 104.0, 106.0, 108.0, 110.0, 112.0,
     &114.0, 116.0, 118.0, 120.0, 122.0, 
     &124.0, 126.0, 128.0, 130.0, 132.0, 134.0,
     &136.0, 138.0, 140.0, 142.0, 144.0, 
     &146.0, 148.0, 150.0, 152.0, 154.0, 156.0,
     &160.0, 170.0, 180.0, 190.0, 200.0,
     &210.0, 220.0, 230.0, 240.0, 251.0, 261.0,
     &271.0, 281.0, 291.0, 301.0, 311.0, 
     &321.0, 331.0, 341.0, 351.0, 361.0, 371.0,
     &381.0, 391.0, 401.0, 411.0, 421.0, 
     &431.0, 441.0, 451.0, 461.0, 471.0, 481.0,
     &491.0, 501.0, 511.0, 521.0, 531.0, 
     &541.0, 551.0, 561.0, 571.0, 581.0, 591.0,
     &601.0, 611.0, 621.0, 631.0, 641.0, 
     &651.0, 661.0, 671.0, 681.0, 691.0, 701.0,
     &711.0, 721.0, 731.0, 741.0, 751.0, 
     &761.0, 771.0, 781.0, 791.0, 801.0, 811.0,
     &821.0, 831.0, 841.0, 851.0, 861.0, 
     &871.0, 881.0, 891.0, 901.0, 911.0, 921.0,
     &931.0, 941.0, 951.0, 961.0, 971.0, 
     &981.0, 991.0, 1001.0, 1011.0, 1021.0,
     &1031.0, 1041.0, 1051.0, 1061.0, 1071.0, 
     &1081.0, 1091.0, 1101.0, 1111.0, 1121.0,
     &1131.0, 1141.0, 1151.0, 1161.0, 1171.0, 
     &1181.0, 1191.0, 1201.0, 1211.0, 1221.0,
     &1231.0, 1241.0, 1251.0, 1261.0, 1271.0, 
     &1281.0, 1291.0, 1301.0, 1311.0, 1321.0,
     &1331.0, 1341.0, 1351.0, 1361.0, 1371.0, 
     &1381.0, 1391.0, 1401.0, 1411.0, 1421.0,
     &1431.0, 1441.0, 1451.0, 1461.0, 1471.0, 
     &1481.0, 1491.0, 1501.0, 1511.0, 1521.0,
     &1531.0, 1541.0, 1551.0, 1561.0, 1571.0, 
     &1581.0, 1591.0, 1601.0, 1611.0, 1621.0,
     &1631.0, 1641.0, 1651.0, 1661.0, 1671.0, 
     &1681.0, 1691.0, 1701.0, 1711.0, 1721.0,
     &1731.0, 1741.0, 1751.0, 1761.0, 1771.0, 
     &1781.0, 1791.0, 1801.0, 1811.0, 1821.0,
     &1831.0, 1841.0, 1851.0, 1861.0, 1871.0, 
     &1881.0, 1891.0, 1901.0, 1911.0, 1921.0,
     &1931.0, 1941.0, 1951.0, 1961.0, 1971.0, 
     &1981.0, 1991.0, 2001.0, 2011.0, 2021.0,
     &2031.0, 2041.0, 2051.0, 2061.0, 2071.0, 
     &2081.0, 2091.0, 2101.0, 2111.0, 2121.0,
     &2131.0, 2141.0, 2151.0, 2161.0, 2171.0, 
     &2181.0, 2191.0, 2201.0, 2211.0, 2221.0,
     &2231.0, 2241.0, 2251.0, 2261.0, 2271.0, 
     &2281.0, 2291.0, 2301.0, 2311.0, 2321.0,
     &2331.0, 2341.0, 2351.0, 2361.0, 2371.0, 
     &2381.0, 2391.0, 2401.0, 2411.0, 2421.0,
     &2431.0, 2441.0, 2451.0, 2461.0, 2471.0, 
     &2481.0, 2491.0, 2501.0 ]

      Fim_y= [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     &0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.016031827022947347, 
     &-0.016883499865276434,-0.01433259282685784 ,-0.01510664711383782,
     &-0.018676300595434698, -0.019877950048475326, 
     &-0.018523870999752022,-0.010999026338845752,
     &-0.0071240510082503814, 1.0510039291942015, 3.6997641774029235, 
     &5.42543775249682     ,7.000776481442942    ,8.561316317459609,     
     &10.178603912697717   , 11.905895049610242  , 13.805349352046996,
     &15.966217814958517   ,18.602981007013575   ,22.70443456451701,     
     &26.833499462470712   , 30.016335467160797  , 32.29520365573807,
     &33.76863593152834    ,34.55777829352152    ,34.769576371379,       
     &34.5027637207783     , 33.83746535826642   , 32.85184284860148,
     &31.607384630212398   ,30.15381673575334    ,28.53750641910118,     
     &26.80012979801171    , 24.980050808058184  , 23.11297249488223,
     &21.226953772763107   ,19.34622388476728    ,17.48305219421484,     
     &15.649790844298346   , 13.854068034810796  , 12.102579237235634,
     &10.402839197687022   ,8.757464999782284    ,7.167031099435345,  
     &5.633073412671622    , 4.156791449657249   , 2.738575313067803,
     &1.377296390242001    ,0.07128212447546677  ,-1.1807790389892374,
     &-2.3785043675635227  , -3.5241209229298205 , 
     &-4.620060415377433   ,-5.667232736020825   ,-6.667235686293791,  
     &-7.620698228078953   , -8.531402598886025  , 
     &-9.399881704578938   ,-10.22813506615573   ,-11.017029317080558,   
     &-11.769044954955906  , -12.48717122922285  , 
     &-13.169844118584415  ,-13.819390840463104  ,-14.437558379234437,  
     &-15.026949587421914  , -15.587523825451711 , 
     &-16.121154606999703  ,-16.628359463232364  ,-17.11078149694818,  
     &-17.56915615015241   , -18.004420445748003 , 
     &-18.418435944675725  ,-18.811620790312986  ,-19.184505992707656,   
     &-19.53842115998419   , -19.874179221339986 , 
     &-20.192651628837258  ,-20.494155438068326  ,-20.77999502558053,  
     &-21.049780647308797  , -21.3051365273696   , 
     &-21.546539612068297  ,-21.774724187530868  ,-21.99054101831002,   
     &-22.19401770692825   , -22.385810947963    , -22.566042511448767,
     &-22.736481848792987  ,-22.895763718731217  ,-23.046602551881275, 
     &-23.18714141359533   , -23.318652655984117 , -23.442755480792368,
     &-23.55806139000317   ,-23.6652735974234    ,-23.76511924023524, 
     &-23.857615651961424  , -23.94291945483077  , -24.021747864677796,
     &-24.09433082982495   ,-24.160254604896245  ,-24.22118796693728,    
     &-24.27598267788764   , -24.32572243295854  , -24.370236750889536,
     &-24.410475822912122  ,-24.445578063666364  ,-24.47555419737684,    
     &-24.502088784448272  , -24.524391759056055 , -24.542230383590002,
     &-24.55644783818969   ,-24.56716637815416   ,-24.574223194735822,   
     &-24.57797712203586   , -24.579122697150655 , -24.57660364029591,
     &-24.571118503243444  ,-24.563083555264665  ,-24.5523231599685,     
     &-24.539713009477207  , -24.523394991823842 , -24.50488327536236,
     &-24.4842322265726    ,-24.46150437649652   ,-24.436511893816565,   
     &-24.40913820062609   , -24.38033572508962  , -24.349061398095408,
     &-24.317317108704668  ,-24.283033766219503  ,-24.24769722165324,    
     &-24.212305156708926  , -24.17388433762043  , -24.133720503934054,
     &-24.092416004227932  ,-24.049273732987466  ,-24.00672419055335,    
     &-23.961359026974563  , -23.914919277665454 , -23.868538016082887,
     &-23.820407487029193  ,-23.771527751455725  ,-23.721920275081526,   
     &-23.67159055883348   , -23.619599934414243 , -23.56695107056136,
     &-23.51413925245856   ,-23.459683965078476  ,-23.40627874784636,    
     &-23.34981951320564   , -23.294497056197347 , -23.238167571783755,
     &-23.18221018590598   ,-23.124748162632024  ,-23.067174300550782,   
     &-23.009665959009553  , -22.951109146753964 , -22.89282241437192,
     &-22.832386643169137  ,-22.7735741078483    ,-22.712853275259,      
     &-22.653591057596064  , -22.59165459228115  , -22.53122853962109,
     &-22.470951191829982  ,-22.40913391215972   ,-22.34752177082142,    
     &-22.286160827082764  , -22.22281442622782  , -22.16221467359835,
     &-22.100081148224326  ,-22.037697273679026  ,-21.97528080185179,    
     &-21.912484336592236  , -21.849040025786135 , -21.788207113023017,
     &-21.724562963802217  ,-21.66150696469393   ,-21.597736954937933,   
     &-21.535211518838494  , -21.471818950252555 , -21.41030095268589,
     &-21.347378256463266  ,-21.282679775684578  ,-21.2194687008142,     
     &-21.155917925351368  , -21.0947223067983   , -21.0306163581056,
     &-20.967089312783813  ,-20.905755506043835  ,-20.84381478347673,    
     &-20.78128184469374   , -20.717705417441863 , -20.657267621650625,
     &-20.594788060767033  ,-20.53201629553023   ,-20.46985981753961,    
     &-20.40712763162376   , -20.345525328912604 , -20.285175455941182,
     &-20.223665705529275  ,-20.16102597194358   ,-20.0987262005944,     
     &-20.039289825025975  , -19.976759048723036 , -19.913894034110047,
     &-19.854733698768474  ,-19.792925962330457  ,-19.733665641823084,   
     &-19.671965896889404  , -19.6130537206621   , -19.55245542549631,
     &-19.492568277276604  ,-19.43384833226297   ,-19.37033164929107,    
     &-19.313474987813766  , -19.254160875197627 , -19.195267069576495,
     &-19.135925245688526  ,-19.07696002893518   ,-19.019904256787893,   
     &-18.960761675126587  , -18.903396519313368 , -18.846027805114137,
     &-18.787232632399004  ,-18.73287746755258   ,-18.674085227557114,   
     &-18.61538788072107   , -18.557567856919714 , -18.500824014501497,
     &-18.4418221015764    ,-18.38589406548895   ,-18.327123744890688]

c--- Calculate the 1-loop form factor

      F_1l = mt2*(two-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
     &  /(two*wmass*sinthw)
c--- Anomalous top allowed
      

c--- Calculate the 2-loop form factor        
        
      sqrts = sqrt(s(1,2))
      call lin_interpolate(Fre_x, Fre_y, size(Fre_x), sqrts, Fre)
      call lin_interpolate(Fim_x, Fim_y, size(Fim_x), sqrts, Fim)
      F_2l = dcmplx(Fre,Fim)

      F_2l = MH2/(two*wmass*sinthw)*dcmplx(Fre,Fim)

c--- Wave function renormalisation constant

      dZh  = -w3*dB0h

c--- Assemble the renormalised 2-loop form factor        
 
      ggHmt_c6(2,2) = (9*(2.d0 + c6)*dZh*MH2)/2*F_1l 
      
      ggHmt_c6(2,2) = ggHmt_c6(2,2) + F_2l
      ggHmt_c6(2,2) = (c6*MH2)/(32.d0*Pi**2*vevsq)*ggHmt_c6(2,2)
      ggHmt_c6(1,1)=ggHmt_c6(2,2)*za(1,2)/zb(1,2)
      ggHmt_c6(2,2)=ggHmt_c6(2,2)*zb(1,2)/za(1,2)
      return
      end

      subroutine anomhzzamp_c6(prop34,prop56,za,zb,
     & H4l_c6_gmunu,H4l_c6_qmuqnu)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'zcouple.f'
      include 'AnomTriLinear.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcouple.f'
      include 'qlfirst.f'
      double complex H4l(2,2),prop34,prop56
      double complex H4l_c6_gmunu(2,2),HZZ_c6_gmunu
      double complex H4l_c6_qmuqnu(2,2),HZZ_c6_qmuqnu
      double precision mb2,mt2,sinthw
      double precision qhsq, q1sq, q2sq, MH2, MZ2, dB0h, dZh, Z
      double precision sqrts 
      double complex qlI2,qlI3,C0mt

c--- SM H4l amplitude 
c--- Possible in the future to modify AnomHZZ
      MH2  = hmass**2
      MZ2  = zmass**2
      sinthw=sqrt(xw)

c--- SM Amplitudes for decay
      H4l(1,1)=za(3,5)*zb(4,6)*l1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,1)=za(4,5)*zb(3,6)*r1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(1,2)=za(3,6)*zb(4,5)*l1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,2)=za(4,6)*zb(3,5)*r1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56

c--- c6 correction to HZZ vertex

      if (qlfirst) then
        qlfirst=.false.
      call qlinit 
      endif
    
      qhsq = s(1,2)
      q1sq = s(3,4)
      q2sq = s(5,6)

      dB0h = (-9 + 2*Sqrt(3.)*Pi)/(9.*MH2)
      dZh  = -w2*dB0h 

c--- SM part of HZZ vertex

      HZZ_c6_gmunu = (9*(2.d0 + c6)*dZh*MH2)/2.d0

      HZZ_c6_gmunu = HZZ_c6_gmunu +
     &  3 + ((3*(MH2 - MZ2 + q1sq)*(q1sq - q2sq) + 
     &  3*(MH2 - MZ2 - q1sq)*qhsq)*qlI2(q1sq,MH2,MZ2,1D0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) + 
     &  ((-3*(q1sq - q2sq)*(MH2 - MZ2 + q2sq) - 
     &  3*(-MH2 + MZ2 + q2sq)*qhsq)*qlI2(q2sq,MH2,MZ2,1D0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) + 
     &  ((-3*(q1sq - q2sq)**2 + 3*(-2*MH2 + 2*MZ2 + q1sq + q2sq)*qhsq)*
     &     qlI2(qhsq,MH2,MH2,1D0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) + 
     &  ((6*(MH2 - 2*MZ2)*(q1sq - q2sq)**2 + 
     &       6*(MH2**2 + MZ2**2 + q1sq*q2sq + 3*MZ2*(q1sq + q2sq) - 
     &          MH2*(2*MZ2 + q1sq + q2sq))*qhsq - 6*MZ2*qhsq**2)*
     &     qlI3(qhsq,q1sq,q2sq,MH2,MH2,MZ2,1D0,0))/
     &     (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq))

      HZZ_c6_gmunu = (c6*MH2)/(32.d0*Pi**2*vevsq)*HZZ_c6_gmunu

c---  BSM amplitudes for decay proportional to metric
      
      H4l_c6_gmunu(1,1)= H4l(1,1)*HZZ_c6_gmunu
      H4l_c6_gmunu(2,1)= H4l(2,1)*HZZ_c6_gmunu
      H4l_c6_gmunu(1,2)= H4l(1,2)*HZZ_c6_gmunu
      H4l_c6_gmunu(2,2)= H4l(2,2)*HZZ_c6_gmunu

      
c--- Non-SM part of HZZ vertex         

      HZZ_c6_qmuqnu = (6*(q1sq + q2sq - qhsq))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) - 
     &  (12*MH2*(1 - Log(MH2)))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) + 
     &  (12*MZ2*(1 - Log(MZ2)))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq)) + 
     &  (6*((q1sq - q2sq)*(q1sq*(5*MH2 - 5*MZ2 + q1sq) + 
     &          (MH2 - MZ2 + 5*q1sq)*q2sq) - 
     &  2*(q1sq*(2*MH2 - 2*MZ2 + q1sq) + (-MH2 + MZ2 - 2*q1sq)*q2sq)*
     &  qhsq + (-MH2 + MZ2 + q1sq)*qhsq**2)*qlI2(q1sq,MH2,MZ2,1d0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq))**2 + 
     &  ((-6*(q1sq - q2sq)*((MH2 - MZ2)*q1sq 
     &  + 5*(MH2 - MZ2 + q1sq)*q2sq + 
     &          q2sq**2) + 12*(-(MZ2*q1sq) + MH2*(q1sq - 2*q2sq) + 
     &          2*(MZ2 + q1sq)*q2sq - q2sq**2)*qhsq + 
     &       6*(-MH2 + MZ2 + q2sq)*qhsq**2)*qlI2(q2sq,MH2,MZ2,1d0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq))**2 + 
     &  ((-6*(q1sq - q2sq)**2*(2*MH2 - 2*MZ2 + q1sq + q2sq) + 
     &       12*(q1sq**2 - 4*q1sq*q2sq + q2sq**2 - MH2*(q1sq + q2sq) + 
     &          MZ2*(q1sq + q2sq))*qhsq - 
     &          6*(-4*MH2 + 4*MZ2 + q1sq + q2sq)*
     &    qhsq**2)*qlI2(qhsq,MH2,MH2,1d0,0))
     &    /(q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq))**2 + 
     &  (12*((q1sq - q2sq)**2*(MH2**2 + (MZ2 - q1sq)*(MZ2 - q2sq) + 
     &          2*MH2*(-MZ2 + q1sq + q2sq)) + 
     &       (q1sq*(MH2**2 + MZ2*(MZ2 + q1sq) - 2*MH2*(MZ2 + 2*q1sq)) + 
     &    ((MH2 - MZ2)**2 + 4*MH2*q1sq - 6*MZ2*q1sq + q1sq**2)*q2sq + 
     &          (-4*MH2 + MZ2 + q1sq)*q2sq**2)*qhsq + 
     &    (-2*MH2**2 - 2*MZ2**2 - 2*q1sq*q2sq + MZ2*(q1sq + q2sq) + 
     &          2*MH2*(2*MZ2 + q1sq + q2sq))*qhsq**2 - MZ2*qhsq**3)*
     &     qlI3(qhsq,q1sq,q2sq,MH2,MH2,MZ2,1d0,0))/
     &   (q1sq**2 + (q2sq - qhsq)**2 - 2*q1sq*(q2sq + qhsq))**2

      HZZ_c6_qmuqnu = (c6*MH2)/(32.d0*Pi**2*vevsq)*HZZ_c6_qmuqnu

c---  BSM amplitudes for decay proportional to momenta. See Appendix A of https://arxiv.org/pdf/1902.04756.pdf
      
      H4l_c6_qmuqnu(1,1)=-(za(3,5)*zb(4,5)+za(3,6)*zb(4,6))
     &        *(za(3,5)*zb(3,6)+za(4,5)*zb(4,6))*l1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56*HZZ_c6_qmuqnu   
      H4l_c6_qmuqnu(2,1)=-(za(4,5)*zb(3,5)+za(4,6)*zb(3,6))
     &        *(za(3,5)*zb(3,6)+za(4,5)*zb(4,6))*r1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56*HZZ_c6_qmuqnu    
      H4l_c6_qmuqnu(1,2)=-(za(3,5)*zb(4,5)+za(3,6)*zb(4,6))
     &        *(za(3,6)*zb(3,5)+za(4,6)*zb(4,5))*l1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56*HZZ_c6_qmuqnu         
      H4l_c6_qmuqnu(2,2)=-(za(4,5)*zb(3,5)+za(4,6)*zb(3,6))
     &        *(za(3,6)*zb(3,5)+za(4,6)*zb(4,5))*r1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56*HZZ_c6_qmuqnu      


      return
      end

      subroutine SMggHmtvertex(za,zb,ggHmt)
      implicit none
      include 'masses.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      double complex ggHmt(2,2),C0mt,qlI3
      double precision mt2,sinthw
      
c--- squared masses and sin(thetaw)     
      mt2=mt**2
      sinthw=sqrt(xw)

c----- loop intgral for top
      
      C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)

c------ top quark in the loop
      ggHmt(2,2)=mt2*(two-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
      ggHmt(1,1)=ggHmt(2,2)*za(1,2)/zb(1,2)
      ggHmt(2,2)=ggHmt(2,2)*zb(1,2)/za(1,2)
      
c------ print *,"Here", ggHmt(1,1)
      return
      end
      