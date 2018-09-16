      subroutine sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
      implicit none
      include 'types.f'
c--- set up the necessary parameters for a Standard Model Higgs boson
c---   hwidth : either the NLO value from Spira,
c---                or the LO value calculated here
c---   br,wwbr,zzbr,tautaubr : the LO calculated values
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'couple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'kprocess.f' 
      include 'verbose.f'
      include 'cpscheme.f'
      include 'msbarmasses.f'
      include 'kpart.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      real(dp):: br,gamgambr,zgambr,wwbr,zzbr,tautaubr,x_w,x_z,
     & msqgamgam,hzgamwidth,msqhbb,msqhtautau,
     & pw_bb,pw_tautau,pw_gamgam,pw_ww,pw_zz,pw_zgam,
     & br_sp,tautaubr_sp,gamgambr_sp,wwbr_sp,zzbr_sp,zgambr_sp,
     & GammaHbb0,GammaHbb1,massfrun
      logical:: spira
      common/spira/spira

******************************* H->bb PARAMETERS *******************************
c---- if true, uses the input b-mass in the Yukawa coupling and normalizes
c---- the H->bb cross-section according to the best value from Spira
c---- if false, instead runs the b-mass to the Higgs mass and uses the LO Br

      if(hdecaymode=='bqba') then 
         FixBrHbb=.true.
      else
         FixBrHbb=.false.
      endif

c---  mb_eff is the bottom mass used in the Yukawa coupling
      if (FixBrHbb) then
c--- fix mb to to value in input file
        mb_eff=mb
      else
c--- run mb to appropriate scale
        if (kpart==klord) then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
        endif
      endif

***************************** COMPUTE PARTIAL WIDTHS ***************************
*                                                                              *
*  Included decays are:                                                        *
*                                                                              *
*    H -> bb, H -> tau tau                                                     *
*    H -> WW, H -> ZZ, H -> gamma gamma                                        *
*                                                                              *
********************************************************************************

c--- Compute partial width of H -> bb : note that the mass is mostly retained
c--- in the coupling only, since mb=0 usually (the mass in the phase space)
      pw_bb=msqhbb(hmass**2)/(16._dp*pi*hmass)
     &     *sqrt(1._dp-4._dp*mb**2/hmass**2)

c--- Compute partial width of H -> tau^- tau^+
      pw_tautau=msqhtautau(hmass**2)/(16._dp*pi*hmass)
     &     *sqrt(1._dp-4._dp*mtau**2/hmass**2)

      x_w=4._dp*wmass**2/hmass**2
      x_z=4._dp*zmass**2/hmass**2
c--- Compute partial width of H -> WW
      if (x_w < 1._dp) then
        pw_ww=gwsq/64._dp/pi*hmass**3/wmass**2
     &       *sqrt(1._dp-x_w)*(1._dp-x_w+0.75_dp*x_w**2)
      else
        pw_ww=0._dp
      endif
c--- Compute partial width of H -> ZZ
      if (x_z < 1._dp) then
        pw_zz=gwsq/128._dp/pi*hmass**3/wmass**2
     &       *sqrt(1._dp-x_z)*(1._dp-x_z+0.75_dp*x_z**2)
      else
        pw_zz=0._dp
      endif

c--- Compute partial width of H -> gamma gamma
      pw_gamgam=msqgamgam(hmass)/(16._dp*pi*hmass)
      
c--- Compute partial width of H -> Z gamma
      pw_zgam=hzgamwidth(hmass)
      
****************************** COMPUTE TOTAL WIDTHS ****************************

c--- Calculate the value of hwidth, depending on "spira"
c--- note that the branching ratios "br_sp","gamgambr_sp","wwbr_sp","zzbr_sp"
c--- are not actually used in our calculations   
      if (spira) then
        call higgsp(br_sp,tautaubr_sp,gamgambr_sp,zgambr_sp,wwbr_sp,
     &              zzbr_sp)
      else
        hwidth=pw_bb+pw_tautau+pw_ww+pw_zz+pw_gamgam+pw_zgam
      endif 

c--- complex pole scheme, if desired
      CPscheme=.false.
      if (CPscheme) then
        call interpolate_hto(hmass,hwidth)
c        call hto_main_cpH(hmass,hwidth)
c        write(6,*) 'Call to HTO currently not implemented'
c        stop
      endif
c      hwidth=4.17116e-3_dp ! HTO width at 126 GeV


c--- Set up anomalous width of the Higgs boson if required
      if (abs(hwidth_ratio-1._dp) < 1.e-6_dp) then
        anom_Higgs=.false.
      else
        if ( (kcase==kHZZ_tb) .or. (kcase==kHZZint)
     &   .or.(kcase==kHZZHpi) .or. (kcase==kggZZ4l)
     &   .or.(kcase==kHZZqgI) .or. (kcase==kHWW_tb)
     &   .or.(kcase==kHWWint) .or. (kcase==kHWWHpi)
     &   .or.(kcase==kggWW4l) .or. (kcase==kHVV_tb)
     &   .or.(kcase==kggVV4l) .or. (kcase==kHZZ_jj)
     &   .or.(kcase==kHZZpjt) .or. (kcase==kqqWWqq)
     &   .or.(kcase==kqqZZqq) .or. (kcase==kqqWWss)
     &   .or.(kcase==kqqWZqq)) then
          anom_Higgs=.true.
          keep_SMhiggs_norm=.true.
          hwidth=hwidth*hwidth_ratio
          chi_higgs=hwidth_ratio**(0.25_dp)
          pw_bb=pw_bb*chi_higgs**2
          pw_tautau=pw_tautau*chi_higgs**2
          pw_ww=pw_ww*chi_higgs**2
          pw_gamgam=pw_gamgam*chi_higgs**2
          pw_zgam=pw_zgam*chi_higgs**2
          pw_zz=pw_zz*chi_higgs**2
        else
          write(6,*) 'Anomalous Higgs width not supported for'
          write(6,*) 'this process.'
          stop
        endif
      endif
      
c---- This code now deprecated 
c!===== read in anom_higgs.DAT 
c      call read_anom_higgsp
c!========= read in anomalous higgs couplings for width studies 
c      if(anom_Higgs) then 
c         chi_higgs=hwidth_ratio**(0.25_dp)
c         hwidth=hwidth*chi_higgs**4 
c         if(keep_smhiggs_norm) then 
c            pw_bb=pw_bb*chi_higgs**2 
c            pw_tautau=pw_tautau*chi_higgs**2 
c            pw_ww=pw_ww*chi_higgs**2 
c            pw_gamgam=pw_gamgam*chi_higgs**2 
c            pw_zgam=pw_zgam*chi_higgs**2 
c            pw_zz=pw_zz*chi_higgs**2 
c         endif
c      endif


c      call higgsp(br_sp,tautaubr,gamgambr_sp,zgambr_sp,wwbr_sp,zzbr_sp)
c      write(6,*) 'sethparam:Spira bb     BR',br_sp
c      write(6,*) 'sethparam:MCFM  bb     BR',pw_bb/hwidth
c      write(6,*) 'sethparam:Spira tautau BR',tautaubr_sp
c      write(6,*) 'sethparam:MCFM  tautau BR',pw_tautau/hwidth
c      write(6,*) 'sethparam:Spira Zgam   BR',zgambr_sp
c      write(6,*) 'sethparam:MCFM  Zgam   BR',pw_zgam/hwidth
c      write(6,*) 'sethparam:Spira gamgam BR',gamgambr_sp
c      write(6,*) 'sethparam:MCFM  gamgam BR',pw_gamgam/hwidth
c      write(6,*) 'sethparam:Spira ZZ     BR',zzbr_sp
c      write(6,*) 'sethparam:MCFM  ZZ     BR',pw_zz/hwidth

**************************** COMPUTE BRANCHING RATIOS **************************

c--- Branching ratio H -> bb
      br=pw_bb/hwidth
c--- Branching ratio H -> tau tau
      tautaubr=pw_tautau/hwidth
c--- Branching ratio H -> WW
      wwbr=pw_ww/hwidth
c--- Branching ratio H -> ZZ
      zzbr=pw_zz/hwidth
c--- Branching ratio H -> gamgam
      gamgambr=pw_gamgam/hwidth
c--- Branching ratio H -> Zgam
      zgambr=pw_zgam/hwidth

******************************* H->bb PARAMETERS *******************************

      if (FixBrHbb) then
        br=br_sp
c--- most accurate theoretical calculation of BR
        GamHbb=br*hwidth
c--- BR at LO
        GamHbb0=GammaHbb0(hmass**2,mb**2)
c--- BR at NLO
        GamHbb1=GammaHbb1(hmass**2,mb**2)
      endif

*************************** WRITE OUT BRANCHING RATIOS *************************


      if (verbose) then
       write(6,99) hmass,hwidth,br,tautaubr,wwbr,zzbr,gamgambr,zgambr
       if (spira) then
       write(6,*) '*                                                  *'
       write(6,*) '* Note: branching ratios reported here can be > 1  *'
       write(6,*) '*       since the total Higgs width is calculated  *'
       write(6,*) '*       at NLO but the BR calculation uses a       *'
       write(6,*) '*       partial width at LO only.                  *'
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       endif
      endif

      return

 99   format(/,
     &       ' ****************** Higgs parameters ****************'/, 
     &       ' *                                                  *'/, 
     &       ' *   mass(H) = ',f7.2,'      width(H) = ',e12.5,' *'/,
     &       ' *                                                  *'/, 
     &       ' *              Br( H -> b bbar)  = ',f9.5,'       *'/,
     &       ' *              Br( H -> tau tau) = ',f9.5,'       *'/,
     &       ' *              Br( H -> W W)     = ',f9.5,'       *'/,
     &       ' *              Br( H -> Z Z)     = ',f9.5,'       *'/,
     &       ' *              Br( H -> gam gam) = ',f9.5,'       *'/,
     &       ' *              Br( H -> Z gam)   = ',f9.5,'       *'/,
     &       ' ****************************************************')

      end
