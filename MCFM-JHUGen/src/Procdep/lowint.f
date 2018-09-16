      function lowint(r,wgt)
      implicit none
      include 'types.f'
      real(dp):: lowint
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'kprocess.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
      include 'stopscales.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'runstring.f'
      include 'energy.f'
      include 'VVstrong.f'
      include 'dm_params.f'
      include 'initialscales.f'
      include 'toploopgaga.f'
c---- SSbegin
      include 'reweight.f'
c---- SSend
c --- DSW. To store flavour information :
      include 'nflav.f'
c --- DSW.
      include 'x1x2.f'
      include 'bypart.f'
      integer:: pflav,pbarflav
c--- To use VEGAS random number sequence :
      real(dp):: ran2
      integer:: ih1,ih2,j,k,nvec,sgnj,sgnk,ii,i1,i2,i3,i4
      integer:: i,t
      real(dp):: r(mxdim),W,xmsq,val,val2,ptmp,
     & fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     & pswt,rscalestart,fscalestart,
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     & fxb1(-nf:nf),fxb2(-nf:nf),xmsq_array(-nf:nf,-nf:nf)
      real(dp):: wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      real(dp):: msq1(-nf:nf,-nf:nf),
     & msq4(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      real(dp):: xmsq_bypart(-1:1,-1:1)
      logical:: bin,includedipole,checkpiDpjk
      real(dp):: b1scale,q2scale,q1scale,b2scale
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec,
     & qqb_Waa,qqb_Waa_mad
     & qqb_Zbbmas,qqb_Zbbmas,qqb_totttZ,qqb_totttZ_mad
      common/density/ih1,ih2
      common/bin/bin
      common/BrnRat/BrnRat
      common/bqscale/b1scale,q2scale,q1scale,b2scale
      external qq_tchan_ztq,qq_tchan_ztq_mad
      external qq_tchan_htq,qq_tchan_htq_mad,qq_tchan_htq_amp
      external qqb_gamgam_g,qqb_gmgmjt_gvec
!$omp threadprivate(/bqscale/)

!$omp atomic
      ntotshot=ntotshot+1
      lowint=0._dp
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp
      
      W=sqrts**2
      p(:,:)=0._dp

      call gen_lops(r,p,pswt,*999)

      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      if (kcase.ne.kvlchk6 .and. kcase.ne.ktautau) then
c        call masscuts(p,*999)
c      endif

c---- reject event if any s(i,j) is too small
c      call smalls(s,npart,*999)
c----reject event if any tau is too small
      call smalltau(p,npart,*999)
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

c      call writeout(p)
c      stop
      if (dynamicscale) call scaleset(initscale,initfacscale,p)
      
      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)

c--- Calculate the required matrix elements      
      if     (kcase==kW_only) then
        call qqb_w(p,msq)
      elseif (kcase==kW_1jet) then
        call qqb_w_g(p,msq)
c        call qqb_w_gbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k) .ne. 0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        stop
      elseif (kcase==kWgamma) then
        call qqb_wgam(p,msq)
      elseif (kcase==kWgajet) then
         call qqb_wgam_g(p,msq)
      elseif (kcase==kWbfrmc) then
        call qqb_wbfromc(p,msq)
      elseif (kcase==kW_cjet) then
        call qqb_w_cjet(p,msq)
      elseif (kcase==kWcjet0) then
        call qqb_w_cjet_massless(p,msq)
      elseif (kcase==kWbbmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbjem) then
        call qqb_wbbm_g(p,msq)
      elseif (kcase==kWttmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbbar) then
        call qqb_wbb(p,msq)
      elseif (kcase==kW_2jet) then
        call qqb_w2jet(p,msq)
      elseif (kcase==kW_3jet) then
        call qqb_w2jet_g(p,msq)
      elseif (kcase==kWbbjet) then
        call qqb_wbb_g(p,msq)
      elseif (kcase==kZ_only) then
c        if (ewcorr) then
c          call qqb_z_ew(p,msq)
c        else
          call qqb_z(p,msq)
c        endif
      elseif (kcase==kZ_1jet) then
        call qqb_z1jet(p,msq)
c        call qqb_z1jetbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k) .ne. 0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        pause
      elseif (kcase==kZ_2jet) then
        call qqb_z2jet(p,msq)
      elseif (kcase==kZ_3jet) then
        call qqb_z2jet_g(p,msq)
      elseif (kcase==kZgamma) then
        call qqb_zgam(p,msq)
      elseif (kcase==kZ_2gam) then
        call qqb_zaa(p,msq)
      elseif (kcase==kW_2gam) then
c        if (checkpiDpjk(p)) goto 999
        call qqb_Waa(p,msq)
      elseif (kcase==kZgajet) then
        call qqb_zaj(p,msq)
c        call qqb_zgam_g(p,msq)      ! Old routine
      elseif (kcase==kZ2gajt) then
        call qqb_zaa_g(p,msq)
      elseif (kcase==kZga2jt) then
        call qqb_zaj_g(p,msq)
      elseif (kcase==kZbbmas) then
        call qqb_zbbm(p,msq)
      elseif (kcase==kZbbbar) then
        call qqb_zbb(p,msq)
      elseif (kcase==kZbbjet) then
        call qqb_zbb_g(p,msq)
      elseif (kcase==kWWqqbr) then
        call qqb_ww(p,msq)
      elseif (kcase==kWWnpol) then
        call qqb_ww_unpol(p,msq)
      elseif (kcase==kWW_jet) then
        call qqb_ww_g(p,msq)
      elseif (kcase==kWpWp2j) then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpWp3j) then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)


      elseif (kcase==kWW_jet) then
        call qqb_ww_g(p,msq)
c--- Check of gvec routine
c       n(1)=1._dp
c       n(2)=0._dp
c       n(3)=0._dp
c       n(4)=0._dp
c       call qqb_ww_gvec(p,n,2,msqn)
c       n(1)=0._dp
c       n(2)=1._dp
c       n(3)=0._dp
c       n(4)=0._dp
c       call qqb_ww_gvec(p,n,2,msqmad)
c--- polarization vectors in general (for p7)
c       n(1)=p(7,2)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=-p(7,1)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=0._dp
c       n(4)=0._dp       
c       call qqb_ww_gvec(p,n,7,msqn)
c       n(1)=p(7,1)*p(7,3)/p(7,4)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=p(7,2)*p(7,3)/p(7,4)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=-sqrt(p(7,1)**2+p(7,2)**2)/p(7,4)
c       n(4)=0._dp       
c       call qqb_ww_gvec(p,n,7,msqmad)
c       do j=-2,2
c       do k=-2,2
c       if (msq(j,k) .ne. 0._dp) write(6,'(2i3,2e18.8,f16.9)') 
c     &    j,k,msq(j,k),msqmad(j,k)+msqn(j,k),
c     &    msq(j,k)/(msqmad(j,k)+msqn(j,k))
c       enddo
c       enddo
c       pause
c--- Madgraph check
c        call qqb_ww_g_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msq(j,k) .ne. 0._dp)
c     &    write(6,'(2i3,2e18.8,f16.9)')
c     &    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause

c      elseif (kcase==kWW2jet) then
c        call qqb_wwg_g(p,msq)

c        call qqb_ww_gg_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msqmad(j,k) .ne. 0._dp) write(6,'(2i3,2e18.8,f16.9)') 
c     &    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause
      elseif (kcase==kWpWp2j) then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpWp3j) then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)
      elseif (kcase==kZZlept) then
        call qqb_zz(p,msq)
      elseif (kcase==kZZ_jet) then
        call qqb_zz_g(p,msq)
      elseif (kcase==kWHbbar) then
        call qqb_wh(p,msq)
      elseif (kcase==kWH1jet) then
        call qqb_WH1jet(p,msq)
      elseif (kcase==ktwojet) then
        call qqb_twojet(p,msq)
c      elseif (kcase==kthrjet) then
c        call qqb_3jet(p,msq)
      elseif (kcase==kdirgam) then
        call qqb_dirgam(p,msq)
      elseif (kcase==khflgam) then
        call qqb_hflgam(p,msq)
      elseif (kcase==kgamgam) then
        call qqb_gamgam(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        call qqb_gamgam_mad(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        pause
      elseif (kcase==kgg2gam) then
         if(toploopgaga) then
            msq(:,:)=zip
!            call gg_2gam(p,msq)
            call gggaga_mt(p,msq(0,0))
         else
            call gg_2gam(p,msq)
         endif
      elseif (kcase==kgmgmjt) then
c      call checkgvec(+2, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0,-1,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0, 2,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 1,-1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
!        call qqb_gamgam_g(p,msq)
         call qqb_gmgmjt(p,msq)
      elseif (kcase==ktrigam) then
        call qqb_trigam(p,msq)
      elseif (kcase==kfourga) then 
         call qqb_fourgam(p,msq)
      elseif (kcase==kgamjet) then
        call qqb_dirgam_g(p,msq)
      elseif (kcase==kWH__WW) then
        call qqb_wh_ww(p,msq)
      elseif (kcase==kWH__ZZ) then
        call qqb_wh_zz(p,msq)
      elseif (kcase==kWHgaga) then
        call qqb_wh_gaga(p,msq) 
      elseif (kcase==kZHbbar) then
        call qqb_zh(p,msq)
      elseif (kcase==kZH1jet) then
        call qqb_ZH1jet(p,msq)
      elseif (kcase==kZH__WW) then
        call qqb_zh_ww(p,msq)
      elseif (kcase==kZH__ZZ) then
        call qqb_zh_zz(p,msq)
      elseif (kcase==kZHgaga) then
        call qqb_zh_gaga(p,msq)
      elseif (kcase==kggfus0) then
        call gg_h(p,msq)
      elseif (kcase==kHigaga) then
        call gg_hgamgam(p,msq)
      elseif (kcase==kHi_Zga) then
        call gg_hzgam(p,msq)
      elseif (kcase==kHWW_4l) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW2lq) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW_tb) then
        call qqb_hww_tb(p,msq)
      elseif ((kcase==kHWWint) .or. (kcase==kHWWHpi)
     &   .or. (kcase==kggWW4l)) then
        call gg_ww_int(p,msq)
      elseif (kcase==kggWWbx) then
        msq(:,:)=zip
        call gg_WW(p,msq(0,0))
      elseif (kcase==kHZZ_4l) then
        call qqb_hzz(p,msq)
      elseif (kcase==kHZZ_tb) then
        call gg_hzz_tb(p,msq)
      elseif (kcase==kHVV_tb) then
        call gg_hvv_tb(p,msq)
      elseif (kcase==kggVV4l) then
        call gg_VV_all(p,msq)
      elseif (kcase==kggVVbx) then
        msq(:,:)=0._dp
        call gg_VV(p,msq(0,0))
      elseif (kcase==kHZZint) then
        call gg_zz_int(p,msq)
      elseif (kcase==kHZZHpi) then
        call gg_zz_Hpi(p,msq)
      elseif (kcase==kggZZ4l) then
        call gg_zz_all(p,msq)
      elseif (kcase==kggZZbx) then
        msq(:,:)=0._dp
        call gg_ZZ(p,msq(0,0))
      elseif (kcase==kHZZqgI) then 
         call qg_Hint_ZZ(p,msq)
      elseif (kcase==kH_1jet) then
        call qqb_hg(p,msq)
      elseif (kcase==kttZbbl) then
        call qqbZtt(p,msq)
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
        call qqb_QQbdk(p,msq)
      elseif (kcase==ktt_bbu) then
        call qqb_QQbdku(p,msq)
      elseif (kcase==kqq_ttg) then
       call qqb_QQbdk_g(p,msq)
      elseif (kcase==ktt_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==kbb_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==kcc_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==ktt_glu) then
         call qqb_QQb_g(p,msq)
      elseif (kcase==kbq_tpq) then
        call bq_tpq(p,msq)
      elseif (kcase==kttdkay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kt_bbar) then
      call qqb_tbbdk(p,msq)
      elseif (kcase==ktdecay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kW_tndk) then
        call qqb_w_tndk(p,msq)
      elseif (kcase==kW_twdk) then
        call qqb_w_twdk(p,msq)
      elseif (kcase==kWtbwdk) then
        call qqb_wtbwdk(p,msq)
      elseif (kcase==kWtbndk) then
        call qqb_wtbndk(p,msq)
      elseif (kcase==ktottth) then
        call qqb_tottth(p,msq)
      elseif (kcase==kqq_tth) then
        call qqb_tth(p,msq)
      elseif (kcase==ktth_ww) then
        call qqb_tth(p,msq)
      elseif (kcase==kqq_ttz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqqtthz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqq_ttw) then
        call qqb_ttw(p,msq)
      elseif (kcase==khttjet) then
        call qqb_higgs(p,msq)
      elseif (kcase==kggfus1) then
        call gg_hg(p,msq)
      elseif (kcase==kHgagaj) then
        call gg_hgagag(p,msq)
      elseif (kcase==kHWWjet) then
        call gg_hWWg(p,msq)
      elseif (kcase==kHZZjet) then
        call gg_hZZg(p,msq)
      elseif (kcase==kHWW2jt) then
        call gg_hWWgg(p,msq)
      elseif (kcase==kHZZ2jt) then
        call gg_hZZgg(p,msq)
      elseif (kcase==kHWW3jt) then
        call gg_hWWggg(p,msq)
      elseif (kcase==kHZZ3jt) then
        call gg_hZZggg(p,msq)
      elseif (kcase==kattjet) then
        call qqb_higgs_odd(p,msq)
      elseif (kcase==kqq_Hqq) then
        call VV_hqq(p,msq)
      elseif (kcase==kqq_Hgg) then
        call VV_Hgaga(p,msq)
      elseif (kcase==kqqHqqg) then
        call VV_hqq_g(p,msq)
      elseif (kcase==kqq_HWW) then
        call VV_HWW(p,msq)
      elseif (kcase==kqq_HZZ) then
        call VV_HZZ(p,msq)
      elseif (kcase==ktautau) then
        call qqb_tautau(p,msq)
      elseif (kcase==kqg_tbq) then
        call qg_tbq(p,msq)
c--- Check of gvec routines
c      call checkgvec(+2,0,2,p,qg_tbq,qg_tbq_gvec)
c      call checkgvec(-1,0,2,p,qg_tbq,qg_tbq_gvec)
      elseif (kcase==kqqZZqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_ZZqqstrong(p,msq)
        else
          call qq_ZZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWWqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWqqstrong(p,msq)
        else
          call qq_WWqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqVVqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
c          call qq_VVqqstrong(p,msq)
          write(6,*) 'Not yet implemented'
          stop
        else
          call qq_VVqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWWss) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWssstrong(p,msq)
        else
          call qq_WWss(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWZqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WZqqstrong(p,msq)
        else
          call qq_WZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqgtbqq) then
        call qg_tbq_g(p,msq)
      elseif (kcase==k4ftwdk) then
        call qg_tbqdk(p,msq)
      elseif (kcase==k4ftjet) then
        call qg_tbqdk_g(p,msq)
      elseif (kcase==kqq_tbg) then
        call qq_tbg(p,msq)
c--- Check of gvec routines
c      call checkgvec(2,-1,5,p,qq_tbg,qq_tbg_gvec)
c      call checkgvec(-1,2,5,p,qq_tbg,qq_tbg_gvec)
      elseif (kcase==kqqtbgg) then
        call qq_tbg_g(p,msq)
      elseif (kcase==kepem3j) then
        call epem3j(p,msq)
      elseif (kcase==kgQ__ZQ) then
        call gQ_zQ(p,msq)
      elseif (kcase==kZccmas) then
        call qqb_zccm(p,msq)
      elseif (kcase==kggfus2) then
        call gg_hgg(p,msq)
      elseif (kcase==kgagajj) then
        call gg_hgg(p,msq)
      elseif (kcase==kggfus3) then
        call gg_hggg(p,msq)
      elseif (kcase==kW_bjet) then
        call qqb_wbjet(p,msq)
      elseif (kcase==kWcjetg) then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (kcase==kZ_bjet) then
        call qqb_zbjet(p,msq)
      elseif (kcase==kZbjetg) then
        call qqb_zbjet_g(p,msq)
      elseif (kcase==kH_tjet) then
        call qq_tchan_htq(p,msq)
      elseif (kcase==kH_tdkj) then
        call qq_tchan_htq_dk(p,msq)
      elseif (kcase==kZ_tjet) then
        call qq_tchan_ztq(p,msq)
      elseif (kcase==kZ_tdkj) then
         call qq_tchan_ztq_dk(p,msq)
      elseif (kcase==kZtdk2j) then
         call qq_tchan_ztqg_dk(p,msq)
      elseif (kcase==kZt2jet) then
        call qq_tchan_ztqg(p,msq)
      elseif (kcase==kHHpair) then
        call gg_HH(p,msq)
      elseif (kcase==kdm_jet) then 
         call qqb_dm_monojet(p,msq)
      elseif (kcase==kdm_gam) then 
         call qqb_dm_monophot(p,msq)
      elseif ( kcase==kdm2jet) then 
         call qqb_dm_monojet_g(p,msq) 
      elseif ( kcase==kdm_gaj) then 
         call qqb_dm_monophot_g(p,msq)
      elseif (kcase==kvlchk2) then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
      elseif (kcase==kvlchk3) then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
      elseif (kcase==kvlchk4) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)
      elseif (kcase==kvlchk5) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp
        fx2(-1)=4._dp
      elseif (kcase==kvlchk6) then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp*xx(1)
        fx2(-1)=4._dp/xx(2)
      elseif (kcase==kvlchk8) then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0._dp
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp/xx(1)
        fx2(-1)=6._dp/xx(2)
      elseif (kcase==kvlchkm) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)      
      elseif (kcase==kvlchm3) then
        taumin=(2._dp*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
      elseif ((kcase==kvlchwt) .or. (kcase==kvlchwn)
     &   .or. (kcase==kvlchwg) .or. (kcase==kvlchwh)) then
        taumin=0.0001_dp
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
      else
        write(6,*) 'Unimplemented process in lowint : kcase=',kcase
        stop 
      endif
      
      
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0._dp
      enddo
      enddo 

      currentPDF=0

c--- do not calculate the flux if we're only checking the volume      
c      if (case(1:4) .ne. 'vlch') then      
      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)
c      endif
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0._dp
!      if (PDFerrors) then
!        call InitPDF(currentPDF)
!      endif

c--- calculate PDF's  
      if (((kcase==kqg_tbq) .or. (kcase==k4ftwdk))
     &     .and. (dynamicscale .eqv. .false.)) then
c--- for single top + b, make sure to use two different scales
         if (PDFerrors) then
!$omp critical(PDFerrors)
            call InitPDF(currentPDF)
            call fdist(ih1,xx(1),facscale_H,fx1_H)
            call fdist(ih2,xx(2),facscale_H,fx2_H)
            call fdist(ih1,xx(1),facscale_L,fx1_L)
            call fdist(ih2,xx(2),facscale_L,fx2_L)
!$omp end critical(PDFerrors)
         else
            call fdist(ih1,xx(1),facscale_H,fx1_H)
            call fdist(ih2,xx(2),facscale_H,fx2_H)
            call fdist(ih1,xx(1),facscale_L,fx1_L)
            call fdist(ih2,xx(2),facscale_L,fx2_L)
         endif
      do j=-nf,nf
        if (j == 0) then   ! heavy quark line has gluon init. state
          fx1(j)=fx1_H(j)
          fx2(j)=fx2_H(j)
        else
          fx1(j)=fx1_L(j)
          fx2(j)=fx2_L(j)
        endif
      enddo
      else!if (case(1:4) .ne. 'vlch') then
        if ((kcase==kbq_tpq) .or. (kcase==kqg_tbq)) then   
c--- single top: allow for different scales on each leg  
c---  (applies only if dynstring = 'DDIS')
          if (dynstring == 'DDIS') then
             if (PDFerrors) then
!$omp critical(PDFerrors)
                call InitPDF(currentPDF)
                call fdist(ih1,xx(1),b1scale,fxb1)
                call fdist(ih2,xx(2),q2scale,fx2)
                call fdist(ih1,xx(1),q1scale,fx1)
                call fdist(ih2,xx(2),b2scale,fxb2)
!$omp end critical(PDFerrors)
             else
                call fdist(ih1,xx(1),b1scale,fxb1)
                call fdist(ih2,xx(2),q2scale,fx2)
                call fdist(ih1,xx(1),q1scale,fx1)
                call fdist(ih2,xx(2),b2scale,fxb2)
             endif
        else          
           if (PDFerrors) then
!$omp critical(PDFerrors)
              call InitPDF(currentPDF)
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
!$omp end critical(PDFerrors)
           else
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
           endif
        endif
        else   
c--- usual case            
           if (PDFerrors) then
!$omp critical(PDFerrors)
              call InitPDF(currentPDF)
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
!$omp end critical(PDFerrors)
           else
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
           endif
        endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav    

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j==0).and.(k==0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j==0) .or. (k==0)) goto 20
      endif

      if (omitgg) then 
      if ((j==0) .and. (k==0)) goto 20
      endif

      if     ((kcase==kbq_tpq) .and. (dynstring == 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (abs(j) == 5) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (abs(k) == 5) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      elseif ((kcase==kqg_tbq) .and. (dynstring == 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (j == 0) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (k == 0) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      else
c--- DEFAULT
        xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      endif

      xmsq=xmsq+xmsqjk
      xmsq_array(j,k)=xmsqjk
      
      if     (j > 0) then
        sgnj=+1
      elseif (j < 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k > 0) then
        sgnk=+1
      elseif (k < 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      if (currentPDF == 0) then
        xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk
      endif
      
 20   continue
      enddo
      enddo

      if (currentPDF == 0) then
        lowint=flux*pswt*xmsq/BrnRat
      endif
            
c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        PDFwgt(currentPDF)=flux*pswt*xmsq/BrnRat*wgt/itmx
!$omp atomic
        PDFxsec(currentPDF)=PDFxsec(currentPDF)
     &     +PDFwgt(currentPDF)
        currentPDF=currentPDF+1
        if (currentPDF <= maxPDFsets) goto 777
      endif    

        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     &        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      val=wgt*flux*pswt/BrnRat
      do j=-1,1
         do k=-1,1
!$omp atomic            
        lord_bypart(j,k)=lord_bypart(j,k)+
     &       val*xmsq_bypart(j,k)
      enddo
      enddo


      
      val=lowint*wgt
      val2=val**2
      if(val.ne.val) then
         write(6,*) 'lowint val = ',val
         write(6,*) 'Discarding point with random variables',r
         lowint=zip
         val=zip
         goto 999
      endif
c---  SSbegin
      lowint = lowint*reweight
c---  SSend
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (abs(val) > wtmax)) then
!$omp critical(MaxWgt)
        wtmax=abs(val)
!$omp end critical(MaxWgt)
      endif


      if (bin) then
         call nplotter(pjet,val,val2,0)
c---  POWHEG-style output if requested
         if (writepwg) then
!$omp critical(pwhgplot)
            call pwhgplotter(p,pjet,val,0)
!$omp end critical(pwhgplot)
         endif
      endif

c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = abs(val)
c--- note that r(ndim+2) is reserved for new_pspace in real, so unused at LO
        if (r(ndim+2) < (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs<wtmax) then
            newwt = 1._dp
          else
            newwt = wtabs/wtmax
          endif
          if (newwt > 1.0_dp) then
            write(6,*) 'WARNING : lowint : event with |weight| > 1.',
     &                 ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*sign(one,val)
c         call nplotter(pjet,newwt,newwt,0)
!$omp critical(LowintWriteLHE)
          call mcfm_writelhe(pjet,xmsq_array,newwt)
!$omp end critical(LowintWriteLHE)
        endif
      endif

      return

 999  continue
      lowint=0._dp
!$omp atomic
      ntotzero=ntotzero+1
      
      return
      end


