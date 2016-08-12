      !-- H -> [l^-(i3) lb(i4)] [l^-(i5) lb(i6)]
      function anomhzzamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'mxpart.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhzzamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4
      double complex ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3

c------ HZZ DECAY CONVENTIONS
      IF( AllowAnomalousCouplings ) THEN

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

      aa1 =ghz1_dyn*zmass**2/shat
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaBSM**2
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
      

      !-- H -> [l^-(i3) lb(i4)] [l^-(i5) lb(i6)]
      function anomhwwamp(i3,i4,i5,i6,jh,shat,q3_q3,q4_q4,za,zb)
      implicit none
      include 'mxpart.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomhwwamp
      integer i3,i4,i5,i6,jh ! jh is Higgs 1, 2
      double precision shat,q3_q3,q4_q4
      double complex ghw1_dyn,ghw2_dyn,ghw3_dyn,ghw4_dyn
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3

c------ HZZ=-HWW DECAY CONVENTIONS
      IF( AllowAnomalousCouplings ) THEN

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

      aa1 =ghw1_dyn*wmass**2/shat
     &     + (shat-q3_q3-q4_q4)/shat*
     &       (ghw2_dyn
     &       +ghw3_dyn*(shat-q3_q3-q4_q4)/4d0/LambdaBSM**2)
      aa2 =-2d0*ghw2_dyn
     &     -ghw3_dyn*(shat-q3_q3-q4_q4)/2d0/LambdaBSM**2
      aa3 =-2d0*ghw4_dyn

      aa1 = aa1 / wmass**2 
      aa2 = aa2 / wmass**2 
      aa3 = aa3 / wmass**2

      anomhwwamp = ( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )

      ELSE

      anomhwwamp = za(i3,i5)*zb(i4,i6)

      ENDIF

      return

      end
      
