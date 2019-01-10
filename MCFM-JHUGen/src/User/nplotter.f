      subroutine nplotter(p,wt,wt2,nd)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of leptons and jets in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c---     nd:  an integer specifying the dipole number of this contribution
c---          (if applicable), otherwise equal to zero
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'nplot.f'
      include 'nproc.f'

c--- APPLgrid - use of grids
c      include 'ptilde.f'
c      include 'APPLinclude.f'
c--- APPLgrid - end


      double precision p(mxpart,4),wt,wt2
      integer switch,nd
      
c--- This routine simply picks out a process-specific plotting routine
c---  (if available) and falls back to the generic routine otherwise.
c---  So far available: W_only, Z_only, Wbbbar, Wbbmas, WpWp2j, WpWp3j
c---  For the convenience of the user who wants to bail out and do their
c---  own plotting we provide the dummy routine userplotter

c---  the index of the plot - stored in the nplot.f common and used for                                                      
c---  keeping track of the plot index across different files (GPS)                                                           
      nextnplot = 1

c---  first allow for user plots
      call userplotter(p,wt,wt2,nd)

c--- switch:  an integer equal to either 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      if (nd.gt.0) then 
         switch=1
      else
         switch=0
      endif
         
c--- Special plotting routine for WW -> leptons
      if((nproc.eq.61).or.(nproc.eq.126).or.(nproc.eq.127)) then 
         call nplotter_VV(p,wt,wt2,switch,0)
         return 
      endif

      if     (case .eq. 'W_only') then
        call nplotter_W_only(p,wt,wt2,switch)
      elseif (case .eq. 'Z_only') then
        call nplotter_Z_only(p,wt,wt2,switch)
      elseif (case .eq. 'W_cjet') then
        call nplotter_Wbbmas(p,wt,wt2,switch)
      elseif (case .eq. 'Wbbbar') then
        call nplotter_Wbbmas(p,wt,wt2,switch)
      elseif ((case .eq. 'Wbbmas') .or. (case .eq. 'W_bjet')) then
        call nplotter_Wbbmas(p,wt,wt2,switch)
      elseif ((case .eq. 'WpWp2j') .or. (case .eq. 'WpWp3j'))then
        call nplotter_WpWp(p,wt,wt2,switch)
      elseif ((case.eq.'W_1jet') .or. (case.eq.'W_2jet') 
     &   .or. (case.eq.'W_3jet')) then
        call nplotter_Wjets(p,wt,wt2,switch)
      elseif ((case.eq.'HWW_4l') .or. (case.eq.'HWW_tb') 
     &   .or. (case.eq.'HWW2lq') .or. (case.eq.'HWWint') 
     &   .or. (case.eq.'HWWH+i') .or. (case.eq.'ggWW4l')
     &   .or. (case.eq.'WWqqbr') .or. (case.eq.'ggWWbx')) then
        call nplotter_VV(p,wt,wt2,switch,0)
      elseif ((case .eq. 'WW_jet') .or. (case .eq. 'WW2jet')) then
        call nplotter_WW_jet(p,wt,wt2,switch)
c--- photon processes also need to know the dipole number
      elseif ((case.eq.'Wgamma') .or. (case.eq.'Zgamma')
     &   .or. (case.eq.'Wgajet') .or. (case.eq.'Zgajet')) then
        call nplotter_Vgamma(p,wt,wt2,switch,nd)
      elseif (case.eq.'gamgam') then
        call nplotter_gamgam(p,wt,wt2,switch,nd)
      elseif (case .eq. 'gmgmjt') then
        call nplotter_gmgmjt(p,wt,wt2,switch)
      elseif (case.eq.'dirgam') then 
        call nplotter_dirgam(p,wt,wt2,switch,nd)
      elseif (case.eq.'trigam') then 
        call nplotter_trigam(p,wt,wt2,switch)
      elseif (case.eq.'W_2gam')  then
        call nplotter_wgamgam(p,wt,wt2,switch,nd)
      elseif (case.eq.'Z_2gam')  then
        call nplotter_zgamgam(p,wt,wt2,switch,nd)
      elseif (case.eq.'Zgajet')  then
        call nplotter_zgamjet(p,wt,wt2,switch,nd) 
      elseif ((case.eq.'tt_bbl') .or. (case.eq.'tt_ldk')
     &   .or. (case.eq.'tt_bbh') .or. (case.eq.'tt_bbu')
     &   .or. (case.eq.'tt_hdk') .or. (case.eq.'tthWdk')
     &   .or. (case.eq.'tt_udk')) then
        call nplotter_ttbar(p,wt,wt2,switch)
      elseif ((case.eq.'4ftwdk') .or. (case.eq.'dk_4ft')) then
        call nplotter_4ftwdk(p,wt,wt2,switch)
      elseif ((case.eq.'t_bbar') .or. (case.eq.'tdecay')) then
        call nplotter_tbbar(p,wt,wt2,switch)
      elseif (case.eq.'qq_ttw') then
        call nplotter_ttw(p,wt,wt2,switch)
      elseif ((case .eq. 'H_tjet') .or. (case .eq. 'Z_tjet')) then
          call nplotter_Ztj(p,wt,wt2,switch)
      elseif ((case .eq. 'H_tdkj') .or. (case .eq. 'Z_tdkj')) then
           call nplotter_Ztjdk(p,wt,wt2,switch)
      elseif (case .eq. 'qqtthz') then
            call nplotter_ttZ(p,wt,wt2,switch)
      elseif ((case .eq. 'dm_jet').or.(case.eq.'dm2jet')) then 
         call nplotter_dm_monj(p,wt,wt2,switch)
      elseif ((case .eq. 'dm_gam').or.(case.eq.'dm_gaj')) then 
         call nplotter_dm_mongam(p,wt,wt2,switch,nd)
      elseif ((case.eq.'qqZZqq').or.(case.eq.'qqWWqq') 
     &   .or. (case.eq.'qqVVqq').or.(case.eq.'qqWWss')
     &   .or. (case.eq.'qqWZqq').or.(case.eq.'WpmZjj')
     &   .or.(case.eq.'qq_ttg')) then 
         call nplotter_qqZZqq(p,wt,wt2,switch)
      elseif ((case .eq. 'HZZ_4l')
     & .or.   (case .eq. 'HZZ_tb')
     & .or.   (case .eq. 'HZZint')
     & .or.   (case .eq. 'HZZH+i')
     & .or.   (case .eq. 'ggZZ4l') 
     & .or.   (case .eq. 'ggZZbx') 
     & .or.   (case .eq. 'HZZqgI') 
     & .or.   (case .eq. 'ZZlept')
     & .or.   (case .eq. 'ggVV4l') 
     & .or.   (case .eq. 'ggVVbx') 
     & .or.   (case .eq. 'HVV_tb')) then 
         call nplotter_ZZlept(p,wt,wt2,switch)
      else
         call nplotter_auto(p,wt,wt2)
c         call nplotter_generic(p,wt,wt2,switch)
      endif
      
c--- APPLgrid - filling applgrid
c      if (creategrid) call fill_grid(p)
c--- APPLgrid - end

      return
      end
      
