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
      data ewscheme  / +3                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639d-5          /   ! G_F
      data aemmz_inp / 0.0078125d0         /   ! alpha_EM(m_Z)=1/128.89
      data xw_inp    / 0.2312d0            /   ! sin^2(theta_W)
      data wmass_inp / 80.419d0            /   ! W mass
      data zmass_inp / 91.188d0           /   ! Z mass
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

      data mtau,mtausq/1.777d0,3.157729d0/
c----   Note: after v5.6, the masses for top, bottom and charm quarks
c----         are set in the input file

c---  Widths: note that the top width is calculated in the program
c---  The W width of 2.1054 is derived using the measured BR of
c---    10.80 +/- 0.09 % (PDG) and the LO partial width calculation
c---    for Mw=80.398 GeV
      data wwidth,zwidth/2.06d0,2.49d0/
      data tauwidth/2.269d-12/
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
     .   /0.975d0,0.222d0,0.000d0,
     .    0.222d0,0.975d0,0.000d0/
      end
************************************************************************


************************************************************************
*     MS-bar quark massses, mQ_MS(mQ_MS)                               *
*     computed using 1-loop running from pole masses of                *
*     mb=4.75, mt=173.5                                                *
************************************************************************
      block data block_msbarmasses
      implicit none
      include 'msbarmasses.f'
      data mc_msbar/1.275d0/
      data mb_msbar/4.37d0/
      data mt_msbar/166d0/
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

