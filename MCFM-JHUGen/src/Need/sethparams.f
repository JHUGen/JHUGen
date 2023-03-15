      subroutine sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- set up the necessary parameters for a Standard Model Higgs boson
c---   hwidth : either the NLO value from Spira,
c---                or the LO value calculated here
c---   br,wwbr,zzbr,tautaubr : the LO calculated values
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'process.f' 
      include 'verbose.f'
      include 'cpscheme.f'
      include 'widthscheme.f'
      include 'spinzerohiggs_anomcoupl.f'
      double precision br,gamgambr,zgambr,wwbr,zzbr,tautaubr,x_w,x_z,
     & msqgamgam,hzgamwidth,msqhbb,msqhtautau,
     & pw_bb,pw_tautau,pw_gamgam,pw_ww,pw_zz,pw_zgam,
     & br_sp,tautaubr_sp,gamgambr_sp,wwbr_sp,zzbr_sp,zgambr_sp
      logical spira
      common/spira/spira

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
      pw_bb=msqhbb(hmass**2)/(16d0*pi*hmass)
     &     *sqrt(1d0-4d0*mb**2/hmass**2)

c--- Compute partial width of H -> tau^- tau^+
      pw_tautau=msqhtautau(hmass**2)/(16d0*pi*hmass)
     &     *sqrt(1d0-4d0*mtau**2/hmass**2)

      x_w=4d0*wmass**2/hmass**2
      x_z=4d0*zmass**2/hmass**2
c--- Compute partial width of H -> WW
      if (x_w .lt. 1d0) then
        pw_ww=gwsq/64d0/pi*hmass**3/wmass**2
     &       *dsqrt(1d0-x_w)*(1d0-x_w+0.75d0*x_w**2)
      else
        pw_ww=0d0
      endif
c--- Compute partial width of H -> ZZ
      if (x_z .lt. 1d0) then
        pw_zz=gwsq/128d0/pi*hmass**3/wmass**2
     &       *dsqrt(1d0-x_z)*(1d0-x_z+0.75d0*x_z**2)
      else
        pw_zz=0d0
      endif

c--- Compute partial width of H -> gamma gamma
      pw_gamgam=msqgamgam(hmass)/(16d0*pi*hmass)
      
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

c--- complex pole scheme, if desired. 
      CPscheme=.false.
      if (CPscheme) then
        call interpolate_hto(hmass,hwidth)
c        call hto_main_cpH(hmass,hwidth)
c        write(6,*) 'Call to HTO currently not implemented'
c        stop
      endif
c      hwidth=4.17116d-3 ! HTO width at 126 GeV

c--- If the complex pole scheme is not desired, then the width scheme should be selected.
c--- The width scheme has the same enumeration as that of JHUGen. 
c--- 1-> running width, 2-> fixed width, 3 (skipped due to CPscheme as a boolean), 4 -> new running width, 5 -> propagator removal
c--- There is a protection against generation if the resonance is < 2 widths from the decay threshold (i.e. 2*M_Z) for widthscheme 4
c--- ignoreWidthSchemeFourRestriction can be set to true if you want to bypass this protection (at your own risk!)
      widthscheme=2 !the default configuration is no CPscheme and a fixed width
      ignoreWidthSchemeFourRestriction=.false.
      if (widthscheme.eq.4) then
            if ( (hmass - 2d0*zmass < 2d0*hwidth) ) then
                  if (ignoreWidthSchemeFourRestriction) then
                        print *,'ignoring width restriction.'
                  else
                        print *, "reconsider using WidthScheme 4" 
                        print *, "or ignore it"
                        stop 1
                  endif
            endif
            if ( (h2mass .ne. -1d0) ) then
                  if ((h2mass - 2d0*zmass).lt.(2d0*h2width)) then
                        if (ignoreWidthSchemeFourRestriction) then
                              print *,'ignoring width restriction.'
                        else
                              print *, "reconsider using WidthScheme 4" 
                              print *, "or ignore it"
                              stop 1
                        endif
                  endif
            endif
      endif
c--- Set up anomalous width of the Higgs boson if required
      if (abs(hwidth_ratio-1d0) .lt. 1d-6) then
        anom_Higgs=.false.
      else
        if ( (case .eq. 'HZZ_tb') .or. (case .eq. 'HZZint')
     &   .or.(case .eq. 'HZZH+i') .or. (case .eq. 'ggZZ4l')
     &   .or.(case .eq. 'HZZqgI') .or. (case .eq. 'HWW_tb')
     &   .or.(case .eq. 'HWWint') .or. (case .eq. 'HWWH+i')
     &   .or.(case .eq. 'ggWW4l') .or. (case .eq. 'HVV_tb')
     &   .or.(case .eq. 'ggVV4l') .or. (case .eq. 'HZZ_jj')
     &   .or.(case .eq. 'HZZ+jt') .or. (case .eq. 'qqWWqq')
     &   .or.(case .eq. 'qqZZqq') .or. (case .eq. 'qqWWss')
     &   .or.(case .eq. 'qqWZqq')) then
          anom_Higgs=.true.
          keep_SMhiggs_norm=.true.
          hwidth=hwidth*hwidth_ratio
          chi_higgs=hwidth_ratio**(0.25d0)
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
c         chi_higgs=hwidth_ratio**(0.25d0)
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
     .       ' ****************** Higgs parameters ****************'/, 
     .       ' *                                                  *'/, 
     .       ' *   mass(H) = ',f7.2,'      width(H) = ',e12.5,' *'/,
     .       ' *                                                  *'/, 
     .       ' *              Br( H -> b bbar)  = ',f9.5,'       *'/,
     .       ' *              Br( H -> tau tau) = ',f9.5,'       *'/,
     .       ' *              Br( H -> W W)     = ',f9.5,'       *'/,
     .       ' *              Br( H -> Z Z)     = ',f9.5,'       *'/,
     .       ' *              Br( H -> gam gam) = ',f9.5,'       *'/,
     .       ' *              Br( H -> Z gam)   = ',f9.5,'       *'/,
     .       ' ****************************************************')

      end
