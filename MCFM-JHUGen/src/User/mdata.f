      block data electroweak_input
************************************************************************
*     Calculational scheme for EW couplings                            *
************************************************************************
c
c     ewscheme=-1  : Old MCFM default
c                    input values = Gf,alpha(m_Z),m_W,m_Z
c                    output values = sin^2(theta_W),mtop
c
c     ewscheme=0   : Old MadEvent default (= AlpGen with iewopt=2)
c                    input values = sin^2(theta_W),alpha(m_Z),m_Z
c                    output values = m_W,Gf.
c
c     ewscheme=1   : New MCFM default, also Madevent default, "G_mu scheme"
c                    = LUSIFER and AlpGen (iewopt=3) defaults
c                    input values = G_F,m_Z,m_W
c                    output values = sin^2(theta_W),alpha(m_Z).
c
c     ewscheme=2   : input  values = G_F,sin^2(theta_W),alpha(m_Z)
c                    output values = m_W,m_Z.
c
c     ewscheme=3   : User choice. All parameters are left as they are
c                    input here. You have to know what you're doing.
c
      implicit none
      include 'ewinput.f'

c     Begin anomalous couplings
      include 'spinzerohiggs_anomcoupl.f'
      include 'AnomZffCouplings.f'

      data AllowAnomalousCouplings / 1 / ! Disables anomalous couplings
      data AllowAnomalousZffCouplings / 0 / ! Disables anomalous Zff Couplings
      data distinguish_HWWcouplings / 0 /
      data AnomalCouplPR,AnomalCouplDK / 1, 1/
      data channeltoggle_stu / 2 /
      data vvhvvtoggle_vbfvh / 2 /

      data LambdaBSM / 1000d0 /
      data Lambda2BSM / 1000d0 /

c     4th generation quark masses
      data mt_4gen / 100000d0 /
      data mb_4gen / 100000d0 /

c     MARKUS: anomalous couplings for the first resonance
c     Hgg and Hff anomalous couplings
      data kappa_top / (1d0,0d0) / ! SM=1,0
      data kappa_tilde_top / (0d0,0d0) / ! SM=0,0
      data kappa_bot / (1d0,0d0) / ! SM=1,0
      data kappa_tilde_bot / (0d0,0d0) / ! SM=0,0
      data ghg2 / (0d0,0d0) / ! SM=0,0
      data ghg3 / (0d0,0d0) / ! SM=0,0
      data ghg4 / (0d0,0d0) / ! SM=0,0

      data kappa_4gen_top / (0d0,0d0) /
      data kappa_tilde_4gen_top / (0d0,0d0) /
      data kappa_4gen_bot / (0d0,0d0) /
      data kappa_tilde_4gen_bot / (0d0,0d0) /
      data ghg2_4gen / (0d0,0d0) /
      data ghg3_4gen / (0d0,0d0) /
      data ghg4_4gen / (0d0,0d0) /

c     HVV anomalous couplings
      data Lambda_Q  / 10000d0 /
      data Lambda_z1 / 10000d0 /
      data Lambda_z2 / 10000d0 /
      data Lambda_z3 / 10000d0 /
      data Lambda_z4 / 10000d0 /

      data cz_q1sq / 0 /
      data cz_q2sq / 0 /
      data cz_q12sq / 0 /
      data Lambda_z11 / 100d0 /
      data Lambda_z21 / 100d0 /
      data Lambda_z31 / 100d0 /
      data Lambda_z41 / 100d0 /
      data Lambda_z12 / 100d0 /
      data Lambda_z22 / 100d0 /
      data Lambda_z32 / 100d0 /
      data Lambda_z42 / 100d0 /
      data Lambda_z10 / 100d0 /
      data Lambda_z20 / 100d0 /
      data Lambda_z30 / 100d0 /
      data Lambda_z40 / 100d0 /

      data ghz1 / (1d0,0d0) /
      data ghz1_prime / (0d0,0d0) /
      data ghz1_prime2 / (0d0,0d0) /
      data ghz1_prime3 / (0d0,0d0) /
      data ghz1_prime4 / (0d0,0d0) /
      data ghz1_prime5 / (0d0,0d0) /
      data ghz1_prime6 / (0d0,0d0) /
      data ghz1_prime7 / (0d0,0d0) /

      data ghz2 / (0d0,0d0) /
      data ghz2_prime / (0d0,0d0) /
      data ghz2_prime2 / (0d0,0d0) /
      data ghz2_prime3 / (0d0,0d0) /
      data ghz2_prime4 / (0d0,0d0) /
      data ghz2_prime5 / (0d0,0d0) /
      data ghz2_prime6 / (0d0,0d0) /
      data ghz2_prime7 / (0d0,0d0) /

      data ghz3 / (0d0,0d0) /
      data ghz3_prime / (0d0,0d0) /
      data ghz3_prime2 / (0d0,0d0) /
      data ghz3_prime3 / (0d0,0d0) /
      data ghz3_prime4 / (0d0,0d0) /
      data ghz3_prime5 / (0d0,0d0) /
      data ghz3_prime6 / (0d0,0d0) /
      data ghz3_prime7 / (0d0,0d0) /

      data ghz4 / (0d0,0d0) /
      data ghz4_prime / (0d0,0d0) /
      data ghz4_prime2 / (0d0,0d0) /
      data ghz4_prime3 / (0d0,0d0) /
      data ghz4_prime4 / (0d0,0d0) /
      data ghz4_prime5 / (0d0,0d0) /
      data ghz4_prime6 / (0d0,0d0) /
      data ghz4_prime7 / (0d0,0d0) /

      data Lambda_zgs1 / 10000d0 /

      data ghzgs1_prime2 / (0d0,0d0) /
      data ghzgs2 / (0d0,0d0) /
      data ghzgs3 / (0d0,0d0) /
      data ghzgs4 / (0d0,0d0) /

      data ghgsgs2 / (0d0,0d0) /
      data ghgsgs3 / (0d0,0d0) /
      data ghgsgs4 / (0d0,0d0) /

      data Lambda_w1 / 10000d0 /
      data Lambda_w2 / 10000d0 /
      data Lambda_w3 / 10000d0 /
      data Lambda_w4 / 10000d0 /

      data cw_q1sq / 0 /
      data cw_q2sq / 0 /
      data cw_q12sq / 0 /
      data Lambda_w11 / 100d0 /
      data Lambda_w21 / 100d0 /
      data Lambda_w31 / 100d0 /
      data Lambda_w41 / 100d0 /
      data Lambda_w12 / 100d0 /
      data Lambda_w22 / 100d0 /
      data Lambda_w32 / 100d0 /
      data Lambda_w42 / 100d0 /
      data Lambda_w10 / 100d0 /
      data Lambda_w20 / 100d0 /
      data Lambda_w30 / 100d0 /
      data Lambda_w40 / 100d0 /

      data ghw1 / (1d0,0d0) /
      data ghw1_prime / (0d0,0d0) /
      data ghw1_prime2 / (0d0,0d0) /
      data ghw1_prime3 / (0d0,0d0) /
      data ghw1_prime4 / (0d0,0d0) /
      data ghw1_prime5 / (0d0,0d0) /
      data ghw1_prime6 / (0d0,0d0) /
      data ghw1_prime7 / (0d0,0d0) /

      data ghw2 / (0d0,0d0) /
      data ghw2_prime / (0d0,0d0) /
      data ghw2_prime2 / (0d0,0d0) /
      data ghw2_prime3 / (0d0,0d0) /
      data ghw2_prime4 / (0d0,0d0) /
      data ghw2_prime5 / (0d0,0d0) /
      data ghw2_prime6 / (0d0,0d0) /
      data ghw2_prime7 / (0d0,0d0) /

      data ghw3 / (0d0,0d0) /
      data ghw3_prime / (0d0,0d0) /
      data ghw3_prime2 / (0d0,0d0) /
      data ghw3_prime3 / (0d0,0d0) /
      data ghw3_prime4 / (0d0,0d0) /
      data ghw3_prime5 / (0d0,0d0) /
      data ghw3_prime6 / (0d0,0d0) /
      data ghw3_prime7 / (0d0,0d0) /

      data ghw4 / (0d0,0d0) /
      data ghw4_prime / (0d0,0d0) /
      data ghw4_prime2 / (0d0,0d0) /
      data ghw4_prime3 / (0d0,0d0) /
      data ghw4_prime4 / (0d0,0d0) /
      data ghw4_prime5 / (0d0,0d0) /
      data ghw4_prime6 / (0d0,0d0) /
      data ghw4_prime7 / (0d0,0d0) /

c     Form Factors in Higgs Propagator
      data Lambda_ff1 / 1000d0 /
      data Lambda_ff2 / 1000d0 /
      
      data n_ff1 / 2 /
      data n_ff2 / 2 /

c     Higgs tri-linear coupling and WFR constants
      data c6 /15d0/ ! SM = 0 
      data t1_c6 /1d0/ ! All ti should = 1 to include all effects
      data t2_c6 /1d0/
      data t3_c6 /1d0/
      data t4_c6 /1d0/
      data t5_c6 /1d0/
      data t6_c6 /1d0/
      data w1_c6 /1d0/ ! All wi should = 1 to include all effects
      data w2_c6 /1d0/
      data w3_c6 /1d0/
      data w4_c6 /1d0/
      data w5_c6 /1d0/

c     MARKUS: anomalous couplings for 2nd resonance
c     Mass and width for second resonance
      data h2mass  /-1d0/ ! -1 disables the resonance
      data h2width /0d0/

c     Hgg and Hff anomalous couplings
      data kappa2_top / (0d0,0d0) /
      data kappa2_tilde_top / (0d0,0d0) /
      data kappa2_bot / (0d0,0d0) /
      data kappa2_tilde_bot / (0d0,0d0) /
      data gh2g2 / (0d0,0d0) /
      data gh2g3 / (0d0,0d0) /
      data gh2g4 / (0d0,0d0) /

      data kappa2_4gen_top / (0d0,0d0) /
      data kappa2_tilde_4gen_top / (0d0,0d0) /
      data kappa2_4gen_bot / (0d0,0d0) /
      data kappa2_tilde_4gen_bot / (0d0,0d0) /
      data gh2g2_4gen / (0d0,0d0) /
      data gh2g3_4gen / (0d0,0d0) /
      data gh2g4_4gen / (0d0,0d0) /

c     HVV anomalous couplings
      data Lambda2_Q  / 10000d0 /
      data Lambda2_z1 / 10000d0 /
      data Lambda2_z2 / 10000d0 /
      data Lambda2_z3 / 10000d0 /
      data Lambda2_z4 / 10000d0 /

      data c2z_q1sq / 0 /
      data c2z_q2sq / 0 /
      data c2z_q12sq / 0 /
      data Lambda2_z11 / 100d0 /
      data Lambda2_z21 / 100d0 /
      data Lambda2_z31 / 100d0 /
      data Lambda2_z41 / 100d0 /
      data Lambda2_z12 / 100d0 /
      data Lambda2_z22 / 100d0 /
      data Lambda2_z32 / 100d0 /
      data Lambda2_z42 / 100d0 /
      data Lambda2_z10 / 100d0 /
      data Lambda2_z20 / 100d0 /
      data Lambda2_z30 / 100d0 /
      data Lambda2_z40 / 100d0 /

      data gh2z1 / (0d0,0d0) /
      data gh2z1_prime / (0d0,0d0) /
      data gh2z1_prime2 / (0d0,0d0) /
      data gh2z1_prime3 / (0d0,0d0) /
      data gh2z1_prime4 / (0d0,0d0) /
      data gh2z1_prime5 / (0d0,0d0) /
      data gh2z1_prime6 / (0d0,0d0) /
      data gh2z1_prime7 / (0d0,0d0) /

      data gh2z2 / (0d0,0d0) /
      data gh2z2_prime / (0d0,0d0) /
      data gh2z2_prime2 / (0d0,0d0) /
      data gh2z2_prime3 / (0d0,0d0) /
      data gh2z2_prime4 / (0d0,0d0) /
      data gh2z2_prime5 / (0d0,0d0) /
      data gh2z2_prime6 / (0d0,0d0) /
      data gh2z2_prime7 / (0d0,0d0) /

      data gh2z3 / (0d0,0d0) /
      data gh2z3_prime / (0d0,0d0) /
      data gh2z3_prime2 / (0d0,0d0) /
      data gh2z3_prime3 / (0d0,0d0) /
      data gh2z3_prime4 / (0d0,0d0) /
      data gh2z3_prime5 / (0d0,0d0) /
      data gh2z3_prime6 / (0d0,0d0) /
      data gh2z3_prime7 / (0d0,0d0) /

      data gh2z4 / (0d0,0d0) /
      data gh2z4_prime / (0d0,0d0) /
      data gh2z4_prime2 / (0d0,0d0) /
      data gh2z4_prime3 / (0d0,0d0) /
      data gh2z4_prime4 / (0d0,0d0) /
      data gh2z4_prime5 / (0d0,0d0) /
      data gh2z4_prime6 / (0d0,0d0) /
      data gh2z4_prime7 / (0d0,0d0) /

      data Lambda2_zgs1 / 10000d0 /

      data gh2zgs1_prime2 / (0d0,0d0) /
      data gh2zgs2 / (0d0,0d0) /
      data gh2zgs3 / (0d0,0d0) /
      data gh2zgs4 / (0d0,0d0) /

      data gh2gsgs2 / (0d0,0d0) /
      data gh2gsgs3 / (0d0,0d0) /
      data gh2gsgs4 / (0d0,0d0) /

      data Lambda2_w1 / 10000d0 /
      data Lambda2_w2 / 10000d0 /
      data Lambda2_w3 / 10000d0 /
      data Lambda2_w4 / 10000d0 /

      data c2w_q1sq / 0 /
      data c2w_q2sq / 0 /
      data c2w_q12sq / 0 /
      data Lambda2_w11 / 100d0 /
      data Lambda2_w21 / 100d0 /
      data Lambda2_w31 / 100d0 /
      data Lambda2_w41 / 100d0 /
      data Lambda2_w12 / 100d0 /
      data Lambda2_w22 / 100d0 /
      data Lambda2_w32 / 100d0 /
      data Lambda2_w42 / 100d0 /
      data Lambda2_w10 / 100d0 /
      data Lambda2_w20 / 100d0 /
      data Lambda2_w30 / 100d0 /
      data Lambda2_w40 / 100d0 /

      data gh2w1 / (1d0,0d0) /
      data gh2w1_prime / (0d0,0d0) /
      data gh2w1_prime2 / (0d0,0d0) /
      data gh2w1_prime3 / (0d0,0d0) /
      data gh2w1_prime4 / (0d0,0d0) /
      data gh2w1_prime5 / (0d0,0d0) /
      data gh2w1_prime6 / (0d0,0d0) /
      data gh2w1_prime7 / (0d0,0d0) /

      data gh2w2 / (0d0,0d0) /
      data gh2w2_prime / (0d0,0d0) /
      data gh2w2_prime2 / (0d0,0d0) /
      data gh2w2_prime3 / (0d0,0d0) /
      data gh2w2_prime4 / (0d0,0d0) /
      data gh2w2_prime5 / (0d0,0d0) /
      data gh2w2_prime6 / (0d0,0d0) /
      data gh2w2_prime7 / (0d0,0d0) /

      data gh2w3 / (0d0,0d0) /
      data gh2w3_prime / (0d0,0d0) /
      data gh2w3_prime2 / (0d0,0d0) /
      data gh2w3_prime3 / (0d0,0d0) /
      data gh2w3_prime4 / (0d0,0d0) /
      data gh2w3_prime5 / (0d0,0d0) /
      data gh2w3_prime6 / (0d0,0d0) /
      data gh2w3_prime7 / (0d0,0d0) /

      data gh2w4 / (0d0,0d0) /
      data gh2w4_prime / (0d0,0d0) /
      data gh2w4_prime2 / (0d0,0d0) /
      data gh2w4_prime3 / (0d0,0d0) /
      data gh2w4_prime4 / (0d0,0d0) /
      data gh2w4_prime5 / (0d0,0d0) /
      data gh2w4_prime6 / (0d0,0d0) /
      data gh2w4_prime7 / (0d0,0d0) /

c     Form Factors in Higgs Propagator
      data Lambda2_ff1 / 1000d0 /
      data Lambda2_ff2 / 1000d0 /
      
      data n2_ff1 / 0 /
      data n2_ff2 / 0 /

c     anomalous couplings for triple and quartic gauge boson coupling with d_1-3 (SM: set all to 1)
      data dV_A / (1d0,0d0) /
      data dP_A / (1d0,0d0) /
      data dM_A / (1d0,0d0) /
      data dV_Z / (1d0,0d0) /
      data dP_Z / (1d0,0d0) /
      data dM_Z / (1d0,0d0) /
      data dZZWpWm / (1d0,0d0) /
      data dZAWpWm / (1d0,0d0) /
      data dAAWpWm / (1d0,0d0) /
c     anomalous couplings for triple and quartic gauge boson coupling with d_4 (SM: set all to 0)
      data dFour_A / (0d0,0d0) /
      data dFour_Z / (0d0,0d0) /

c     right handed Z couplings to charged leptons (set separately for each Z)
      data reZ / 0d0 / ! SM = (-2*-1*xw_inp)/(2*sqrt(xw_inp*(1-xw_inp)))

c     left handed Z couplings to charged leptons (set separately for each Z)
      data leZ / 0d0 / ! SM = (-1-2*-1*xw_inp)/(2*sqrt(xw_inp*(1-xw_inp)))

c     End anomalous couplings
c     left handed Z couplings to neutrinos 
      data lnZ / 0d0 / ! SM = (1)/(2*sqrt(xw_inp*(1-xw_inp))

c     right handed Z couplings to neutrinos
      data rnZ / 0d0 / ! SM = 0 

c     Anomalous Couplings for the Z to up type quarks in production 
c     Note, these are called as shifts to SM values
      data clanou / 0d0 / ! SM = 0
      data cranou / 0d0 / ! SM = 0

c     Anomalous Couplings for the Z to down type quarks in production
c     Note, these are called as shifts to SM values
      data clanod / 0d0 / ! SM = 0
      data cranod / 0d0 / ! SM = 0

c     Anomalous Couplings for Z to up type quarks in decay 
      data lquZ / 0d0 / ! SM = (tau(j)-two*Q(j)*xw)/two*sqrt(xw*(1d0-xw))
      data rquZ / 0d0 / ! SM = (-two*Q(j)*xw)/two*sqrt(xw*(1d0-xw))

c     Anomalous Couplings for Z to down type quarks in decay
      data lqdZ / 0d0 / ! SM = (tau(j)-two*Q(j)*xw)/two*sqrt(xw*(1d0-xw))
      data rqdZ / 0d0 / ! SM = (-two*Q(j)*xw)/two*sqrt(xw*(1d0-xw))

c     End anomalous couplings


      data ewscheme  / +3                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639d-5          /   ! G_F
      data aemmz_inp / 7.8125d-03 /   ! alpha_EM(m_Z)=1/128
      data xw_inp    / 0.23119d0            /   ! sin^2(theta_W)
      data wmass_inp / 80.399d0            /   ! W mass
      data zmass_inp / 91.1876d0           /   ! Z mass

      end
************************************************************************


      block data transversedefn
************************************************************************
*     Definition to use for computing transverse quantities            *
*         useEt=.false.    transverse momentum [previous default]      *
*         useEt=.true.     transverse energy                           *
************************************************************************
      implicit none
      include 'useet.f'
      data useEt/.false./
      end
************************************************************************


************************************************************************
*     Masses, widths and initial-state flavour information             *
************************************************************************
      block data block_properties
      implicit none
      include 'masses.f'
      include 'nflav.f'
      include 'nores.f'
c--- if true, nores removes all of the gg contribution
      data nores/.false./
c--- Masses: note that "mtausq" is typically used throughout the
c--- program to calculate couplings that depend on the mass, while
c--- "mtau" is the mass that appears in the rest of the matrix
c--- elements and phase space (and may be set to zero in the program,
c--- depending on the process number)

      data mtau,mtausq/1.7768d0,3.15701824d0/
c----   Note: after v5.6, the masses for top, bottom and charm quarks
c----         are set in the input file

c---  Widths: note that the top width is calculated in the program
c---  The W width of 2.1054 is derived using the measured BR of
c---    10.80 +/- 0.09 % (PDG) and the LO partial width calculation
c---    for Mw=80.399 GeV
      data wwidth,zwidth/2.085d0,2.4952d0/
      data tauwidth/2.267d-12/
c--- Number of active flavours in the initial state: this parameter
c--- may be changed in the program for some processes
      data nflav/5/
c--- Masses below here are currently unused
      data md,mu,ms/5d-3,5d-3,1d-1/
      data mel,mmu/0.510997d-3,0.105658389d0/
      end
************************************************************************


************************************************************************
*     CKM matrix entries                                               *
************************************************************************
      block data block_ckm
      implicit none
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      data  Vud  ,  Vus  ,  Vub  ,
     .      Vcd  ,  Vcs  ,  Vcb
     .   /0.974285d0,0.22529d0,0.00362693741331522731d0,
     .    0.225182d0,0.9734244d0,0.0416893794705567033/
c      data  Vud  ,  Vus  ,  Vub  ,
c     .      Vcd  ,  Vcs  ,  Vcb
c     .   /0.97425d0,0.2253d0,0.0413d0,
c     .    0.225d0,0.986d0,0.0411d0/
c      data  Vud  ,  Vus  ,  Vub  ,
c     .      Vcd  ,  Vcs  ,  Vcb
c     .   /1d0,0d0,0.000d0,
c     .    0d0,1d0,0.000d0/

      end
************************************************************************


************************************************************************
*     MS-bar quark massses, mQ_MS(mQ_MS)                               *
*     computed using 1-loop running from pole masses of                *
*     mb=4.66, mt=173.2                                                *
************************************************************************
      block data block_msbarmasses
      implicit none
      include 'msbarmasses.f'
      data mc_msbar/1.275d0/
      data mb_msbar/4.18d0/
      data mt_msbar/160d0/
      end
************************************************************************


************************************************************************
*     Relevant for the H+b process only :                              *
*       susycoup: the deviation of the Higgs coupling from the         *
*                 Standard Model value (S.M. = 1d0)                    *
************************************************************************
      block data block_bH
      implicit none
      include 'susycoup.f'
      data susycoup/1d0/
      end
************************************************************************


************************************************************************
*     Dim. Reg. parameter epsilon, used for checking the proper        *
*      operation of the NLO code in the program                        *
************************************************************************
      block data block_epinv
      implicit none
      include 'epinv.f'
      include 'epinv2.f'
      data epinv/ 1d3/
      data epinv2/1d3/
      end
************************************************************************

