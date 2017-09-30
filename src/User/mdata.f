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
      data ewscheme  / +1                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639d-5          /   ! G_F
      data aemmz_inp / 7.7585538055706d-03 /   ! alpha_EM(m_Z)=1/128.89
      data xw_inp    / 0.2312d0            /   ! sin^2(theta_W)
      data wmass_inp / 80.398d0            /   ! W mass
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
*
*
*
************************************************************************
*     Masses, widths and initial-state flavour information             *
************************************************************************
*
* Moved to mcfm_init.f
*
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
c      data  Vud  ,  Vus  ,  Vub  ,
c     .      Vcd  ,  Vcs  ,  Vcb
c     .   /1d0,0d0,0.000d0,
c     .    0d0,1d0,0.000d0/

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
*
* moved to mcfm_init.f
*


