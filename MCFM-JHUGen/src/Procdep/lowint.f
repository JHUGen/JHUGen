      double precision function lowint(r,wgt)
      use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
      implicit none
      include 'constants.f'
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
      include 'process.f'
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
c---- SSbegin
      include 'reweight.f'
c---- SSend
c --- DSW. To store flavour information :
      include 'nflav.f'
c --- DSW.
      include 'x1x2.f'
      include 'bypart.f'
      integer pflav,pbarflav
c--- To use VEGAS random number sequence :
      double precision ran2
      integer ih1,ih2,j,k,nvec,sgnj,sgnk,ii,i1,i2,i3,i4
      integer i,t
      double precision r(mxdim),W,xmsq,val,val2,ptmp,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt,rscalestart,fscalestart,
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     . fxb1(-nf:nf),fxb2(-nf:nf),xmsq_array(-nf:nf,-nf:nf)
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      double precision msq1(-nf:nf,-nf:nf),
     & msq4(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      double precision xmsq_bypart(-1:1,-1:1)
      logical bin,first,includedipole,checkpiDpjk
      double precision b1scale,q2scale,q1scale,b2scale
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec,
     & qqb_Waa,qqb_Waa_mad
     & qqb_Zbbmas,qqb_Zbbmas,qqb_totttZ,qqb_totttZ_mad
      common/density/ih1,ih2
      common/bin/bin
      common/BrnRat/BrnRat
      common/bqscale/b1scale,q2scale,q1scale,b2scale
      data first/.true./
      save first,rscalestart,fscalestart
      external qq_tchan_ztq,qq_tchan_ztq_mad
      external qq_tchan_htq,qq_tchan_htq_mad,qq_tchan_htq_amp
      external qqb_gamgam_g,qqb_gmgmjt_gvec
!$omp threadprivate(first,rscalestart,fscalestart,/bqscale/)

      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif

      ntotshot=ntotshot+1
      lowint=0d0
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0d0
      
      W=sqrts**2
      p(:,:)=0d0

      call gen_lops(r,p,pswt,*999)

      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      if (case .ne. 'vlchk6' .and. case .ne. 'tautau') then
c        call masscuts(p,*999)
c      endif

c---- reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

c      call writeout(p)
c      stop
      if (dynamicscale) call scaleset(rscalestart,fscalestart,p)
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

c--- Calculate the required matrix elements      
      if     (case .eq. 'W_only') then
        call qqb_w(p,msq)
      elseif (case .eq. 'W_1jet') then
        call qqb_w_g(p,msq)
      elseif (case .eq. 'Wgamma') then
        call qqb_wgam(p,msq)
      elseif (case .eq. 'Wgajet') then
         call qqb_wgam_g(p,msq)
      elseif (case .eq. 'Wbfrmc') then
        call qqb_wbfromc(p,msq)
      elseif (case .eq. 'W_cjet') then
        call qqb_w_cjet(p,msq)
      elseif (case .eq. 'Wcjet0') then
        call qqb_w_cjet_massless(p,msq)
      elseif (case .eq. 'Wbbmas') then
        call qqb_wbbm(p,msq)
      elseif (case .eq. 'Wbbjem') then
        call qqb_wbbm_g(p,msq)
      elseif (case .eq. 'Wttmas') then
        call qqb_wbbm(p,msq)
      elseif (case .eq. 'Wbbbar') then
        call qqb_wbb(p,msq)
      elseif (case .eq. 'W_2jet') then
        call qqb_w2jet(p,msq)
      elseif (case .eq. 'W_3jet') then
        call qqb_w2jet_g(p,msq)
      elseif (case .eq. 'Wbbjet') then
        call qqb_wbb_g(p,msq)
      elseif (case .eq. 'Z_only') then
c        if (ewcorr) then
c          call qqb_z_ew(p,msq)
c        else
          call qqb_z(p,msq)
c        endif
      elseif (case .eq. 'Z_1jet') then
        call qqb_z1jet(p,msq)
      elseif (case .eq. 'Z_2jet') then
        call qqb_z2jet(p,msq)
      elseif (case .eq. 'Z_3jet') then
        call qqb_z2jet_g(p,msq)
      elseif (case .eq. 'Zgamma') then
        call qqb_zgam(p,msq)
      elseif (case .eq. 'Z_2gam') then
        call qqb_zaa(p,msq)
      elseif (case .eq. 'W_2gam') then
c        if (checkpiDpjk(p)) goto 999
        call qqb_Waa(p,msq)
      elseif (case .eq. 'Zgajet') then
        call qqb_zaj(p,msq)
c        call qqb_zgam_g(p,msq)      ! Old routine
      elseif (case .eq. 'Z2gajt') then
        call qqb_zaa_g(p,msq)
      elseif (case .eq. 'Zga2jt') then
        call qqb_zaj_g(p,msq)
      elseif (case .eq. 'Zbbmas') then
        call qqb_zbbm(p,msq)
      elseif (case .eq. 'Zbbbar') then
        call qqb_zbb(p,msq)
      elseif (case .eq. 'Zbbjet') then
        call qqb_zbb_g(p,msq)
      elseif (case .eq. 'WWqqbr') then
        call qqb_ww(p,msq)
      elseif (case .eq. 'WWnpol') then
        call qqb_ww_unpol(p,msq)
      elseif (case .eq. 'WW_jet') then
        call qqb_ww_g(p,msq)
      elseif (case .eq. 'WpWp2j') then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WpWp3j') then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WpmZjj') then
        call qqb_WZjj(p,msq)
      elseif (case .eq. 'WpmZbj') then
        call qqb_WZbj(p,msq)
      elseif (case .eq. 'WpmZbb') then
        call qqb_WZbb(p,msq)
      elseif (case .eq. 'WZbbar') then
        call qqb_wz(p,msq)


      elseif (case .eq. 'WW_jet') then
        call qqb_ww_g(p,msq)
c--- Check of gvec routine
c       n(1)=1d0
c       n(2)=0d0
c       n(3)=0d0
c       n(4)=0d0
c       call qqb_ww_gvec(p,n,2,msqn)
c       n(1)=0d0
c       n(2)=1d0
c       n(3)=0d0
c       n(4)=0d0
c       call qqb_ww_gvec(p,n,2,msqmad)
c--- polarization vectors in general (for p7)
c       n(1)=p(7,2)/dsqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=-p(7,1)/dsqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=0d0
c       n(4)=0d0       
c       call qqb_ww_gvec(p,n,7,msqn)
c       n(1)=p(7,1)*p(7,3)/p(7,4)/dsqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=p(7,2)*p(7,3)/p(7,4)/dsqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=-dsqrt(p(7,1)**2+p(7,2)**2)/p(7,4)
c       n(4)=0d0       
c       call qqb_ww_gvec(p,n,7,msqmad)
c       do j=-2,2
c       do k=-2,2
c       if (msq(j,k) .ne. 0d0) write(6,'(2i3,2e18.8,f16.9)') 
c     .    j,k,msq(j,k),msqmad(j,k)+msqn(j,k),
c     .    msq(j,k)/(msqmad(j,k)+msqn(j,k))
c       enddo
c       enddo
c       pause
c--- Madgraph check
c        call qqb_ww_g_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msq(j,k) .ne. 0d0)
c     .    write(6,'(2i3,2e18.8,f16.9)')
c     .    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause

c      elseif (case .eq. 'WW2jet') then
c        call qqb_wwg_g(p,msq)

c        call qqb_ww_gg_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msqmad(j,k) .ne. 0d0) write(6,'(2i3,2e18.8,f16.9)') 
c     .    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause
      elseif (case .eq. 'WpWp2j') then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WpWp3j') then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (case .eq. 'WpmZjj') then
        call qqb_WZjj(p,msq)
      elseif (case .eq. 'WpmZbj') then
        call qqb_WZbj(p,msq)
      elseif (case .eq. 'WpmZbb') then
        call qqb_WZbb(p,msq)
      elseif (case .eq. 'WZbbar') then
        call qqb_wz(p,msq)
      elseif (case .eq. 'ZZlept') then
        call qqb_zz(p,msq)
      elseif (case .eq. 'ZZ_jet') then
        call qqb_zz_g(p,msq)
      elseif (case .eq. 'WHbbar') then
        call qqb_wh(p,msq)
c      elseif (case .eq. 'twojet') then
c        call qqb_2jet(p,msq)
c      elseif (case .eq. 'thrjet') then
c        call qqb_3jet(p,msq)
      elseif (case .eq. 'dirgam') then
        call qqb_dirgam(p,msq)
      elseif (case .eq. 'hflgam') then
        call qqb_hflgam(p,msq)
      elseif (case .eq. 'gamgam') then
        call qqb_gamgam(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        call qqb_gamgam_mad(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        pause
      elseif (case .eq. 'gmgmjt') then
c      call checkgvec(+2, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0,-1,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0, 2,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 1,-1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
!        call qqb_gamgam_g(p,msq)
         call qqb_gmgmjt(p,msq)
      elseif (case .eq. 'trigam') then
        call qqb_trigam(p,msq)
      elseif (case .eq. 'fourga') then 
         call qqb_fourgam(p,msq)
      elseif (case .eq. 'gamjet') then
        call qqb_dirgam_g(p,msq)
      elseif (case .eq. 'WH__WW') then
        call qqb_wh_ww(p,msq)
      elseif (case .eq. 'WH__ZZ') then
        call qqb_wh_zz(p,msq)
      elseif (case .eq. 'WHgaga') then
        call qqb_wh_gaga(p,msq)
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh(p,msq)
      elseif (case .eq. 'ZH__WW') then
        call qqb_zh_ww(p,msq)
      elseif (case .eq. 'ZH__ZZ') then
        call qqb_zh_zz(p,msq)
      elseif (case .eq. 'ZHgaga') then
        call qqb_zh_gaga(p,msq)
      elseif (case .eq. 'ggfus0') then
        call gg_h(p,msq)
      elseif (case .eq. 'Higaga') then
        call gg_hgamgam(p,msq)
      elseif (case .eq. 'Hi_Zga') then
        call gg_hzgam(p,msq)
      elseif (case .eq. 'HWW_4l') then
        call qqb_hww(p,msq)
      elseif (case .eq. 'HWW2lq') then
        call qqb_hww(p,msq)
      elseif (case .eq. 'HWW_tb') then
        call qqb_hww_tb(p,msq)
      elseif ((case .eq. 'HWWint') .or. (case .eq. 'HWWH+i')
     &   .or. (case .eq. 'ggWW4l')) then
        call gg_ww_int(p,msq)
      elseif (case .eq. 'ggWWbx') then
        msq(:,:)=zip
        call gg_WW(p,msq(0,0))
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz(p,msq)
      elseif (case .eq. 'HZZ_tb') then
        call gg_hzz_tb(p,msq)
      elseif (case .eq. 'HVV_tb') then
        call gg_hvv_tb(p,msq)
      elseif (case .eq. 'ggVV4l') then
        call gg_VV_all(p,msq)
      elseif (case .eq. 'ggVVbx') then
        msq(:,:)=0d0
        call gg_VV(p,msq(0,0))
      elseif (case .eq. 'HZZint') then
        call gg_zz_int(p,msq)
      elseif (case .eq. 'HZZH+i') then
        call gg_zz_Hpi(p,msq)
      elseif (case .eq. 'ggZZ4l') then
        call gg_zz_all(p,msq)
      elseif (case .eq. 'ggZZbx') then
        msq(:,:)=0d0
        call gg_ZZ(p,msq(0,0))
      elseif (case .eq. 'HZZqgI') then 
         call qg_Hint_ZZ(p,msq)
      elseif (case .eq. 'H_1jet') then
        call qqb_hg(p,msq)
      elseif (case .eq. 'ttZbbl') then
        call qqbZtt(p,msq)
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
        call qqb_QQbdk(p,msq)
      elseif (case .eq. 'tt_bbu') then
        call qqb_QQbdku(p,msq)
      elseif (case .eq. 'qq_ttg') then
       call qqb_QQbdk_g(p,msq)
      elseif (case .eq. 'tt_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'bb_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'cc_tot') then
        call qqb_QQb(p,msq)
      elseif (case .eq. 'tt_glu') then
         call qqb_QQb_g(p,msq)
      elseif (case .eq. 'bq_tpq') then
        call bq_tpq(p,msq)
      elseif (case .eq. 'ttdkay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 't_bbar') then
      call qqb_tbbdk(p,msq)
      elseif (case .eq. 'tdecay') then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (case .eq. 'W_tndk') then
        call qqb_w_tndk(p,msq)
      elseif (case .eq. 'W_twdk') then
        call qqb_w_twdk(p,msq)
      elseif (case .eq. 'Wtbwdk') then
        call qqb_wtbwdk(p,msq)
      elseif (case .eq. 'Wtbndk') then
        call qqb_wtbndk(p,msq)
      elseif (case .eq. 'tottth') then
        call qqb_tottth(p,msq)
      elseif (case .eq. 'qq_tth') then
        call qqb_tth(p,msq)
      elseif (case .eq. 'tth_ww') then
        call qqb_tth(p,msq)
      elseif (case .eq. 'qq_ttz') then
        call qqb_ttz(p,msq)
      elseif (case .eq. 'qqtthz') then
        call qqb_ttz(p,msq)
      elseif (case .eq. 'qq_ttw') then
        call qqb_ttw(p,msq)
      elseif (case .eq. 'httjet') then
        call qqb_higgs(p,msq)
      elseif (case .eq. 'ggfus1') then
        call gg_hg(p,msq)
      elseif (case .eq. 'Hgagaj') then
        call gg_hgagag(p,msq)
      elseif (case .eq. 'HWWjet') then
        call gg_hWWg(p,msq)
      elseif (case .eq. 'HZZjet') then
        call gg_hZZg(p,msq)
      elseif (case .eq. 'HWW2jt') then
        call gg_hWWgg(p,msq)
      elseif (case .eq. 'HZZ2jt') then
        call gg_hZZgg(p,msq)
      elseif (case .eq. 'HWW3jt') then
        call gg_hWWggg(p,msq)
      elseif (case .eq. 'HZZ3jt') then
        call gg_hZZggg(p,msq)
      elseif (case .eq. 'attjet') then
        call qqb_higgs_odd(p,msq)
      elseif (case .eq. 'qq_Hqq') then
        call VV_hqq(p,msq)
      elseif (case .eq. 'qq_Hgg') then
        call VV_Hgaga(p,msq)
      elseif (case .eq. 'qqHqqg') then
        call VV_hqq_g(p,msq)
      elseif (case .eq. 'qq_HWW') then
        call VV_HWW(p,msq)
      elseif (case .eq. 'qq_HZZ') then
        call VV_HZZ(p,msq)
      elseif (case .eq. 'tautau') then
        call qqb_tautau(p,msq)
      elseif (case .eq. 'qg_tbq') then
        call qg_tbq(p,msq)
c--- Check of gvec routines
c      call checkgvec(+2,0,2,p,qg_tbq,qg_tbq_gvec)
c      call checkgvec(-1,0,2,p,qg_tbq,qg_tbq_gvec)
      elseif (case .eq. 'qqZZqq') then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_ZZqqstrong(p,msq)
        else
          call qq_ZZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (case .eq. 'qqWWqq') then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWqqstrong(p,msq)
        else
          call qq_WWqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (case .eq. 'qqVVqq') then
c        call getvbfpoint(p)
        if (VVstrong) then
c          call qq_VVqqstrong(p,msq)
          write(6,*) 'Not yet implemented'
          stop
        else
          call qq_VVqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (case .eq. 'qqWWss') then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWssstrong(p,msq)
        else
          call qq_WWss(p,msq)
        endif
c        call comparevbf(msq)
      elseif (case .eq. 'qqWZqq') then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WZqqstrong(p,msq)
        else
          call qq_WZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (case .eq. 'qgtbqq') then
        call qg_tbq_g(p,msq)
      elseif (case .eq. '4ftwdk') then
        call qg_tbqdk(p,msq)
      elseif (case .eq. '4ftjet') then
        call qg_tbqdk_g(p,msq)
      elseif (case .eq. 'qq_tbg') then
        call qq_tbg(p,msq)
c--- Check of gvec routines
c      call checkgvec(2,-1,5,p,qq_tbg,qq_tbg_gvec)
c      call checkgvec(-1,2,5,p,qq_tbg,qq_tbg_gvec)
      elseif (case .eq. 'qqtbgg') then
        call qq_tbg_g(p,msq)
      elseif (case .eq. 'epem3j') then
        call epem3j(p,msq)
      elseif (case .eq. 'gQ__ZQ') then
        call gQ_zQ(p,msq)
      elseif (case .eq. 'Zccmas') then
        call qqb_zccm(p,msq)
      elseif (case .eq. 'ggfus2') then
        call gg_hgg(p,msq)
      elseif (case .eq. 'gagajj') then
        call gg_hgg(p,msq)
      elseif (case .eq. 'ggfus3') then
        call gg_hggg(p,msq)
      elseif (case .eq. 'W_bjet') then
        call qqb_wbjet(p,msq)
      elseif (case .eq. 'Wcjetg') then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (case .eq. 'Z_bjet') then
        call qqb_zbjet(p,msq)
      elseif (case .eq. 'Zbjetg') then
        call qqb_zbjet_g(p,msq)
      elseif (case .eq. 'H_tjet') then
        call qq_tchan_htq(p,msq)
      elseif (case .eq. 'H_tdkj') then
        call qq_tchan_htq_dk(p,msq)
      elseif (case .eq. 'Z_tjet') then
        call qq_tchan_ztq(p,msq)
      elseif (case .eq. 'Z_tdkj') then
         call qq_tchan_ztq_dk(p,msq)
      elseif (case .eq. 'Ztdk2j') then
         call qq_tchan_ztqg_dk(p,msq)
      elseif (case .eq. 'Zt2jet') then
        call qq_tchan_ztqg(p,msq)
      elseif (case .eq. 'HHpair') then
        call gg_HH(p,msq)
      elseif (case .eq. 'dm_jet') then 
         call qqb_dm_monojet(p,msq)
      elseif (case .eq. 'dm_gam') then 
         call qqb_dm_monophot(p,msq)
      elseif ( case.eq.'dm2jet') then 
         call qqb_dm_monojet_g(p,msq) 
      elseif ( case.eq.'dm_gaj') then 
         call qqb_dm_monophot_g(p,msq)
      elseif (case .eq. 'vlchk2') then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      elseif (case .eq. 'vlchk3') then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif (case .eq. 'vlchk4') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)
      elseif (case .eq. 'vlchk5') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0
        fx2(-1)=4d0
      elseif (case .eq. 'vlchk6') then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0*xx(1)
        fx2(-1)=4d0/xx(2)
      elseif (case .eq. 'vlchk8') then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0d0
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=6d0/xx(1)
        fx2(-1)=6d0/xx(2)
      elseif (case .eq. 'vlchkm') then
        taumin=0.0001d0
        bbsqmax=W
        bbsqmin=0d0
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=4d0*xx(1)
        fx2(-1)=2d0/xx(2)      
      elseif (case .eq. 'vlchm3') then
        taumin=(2d0*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=2d0
        fx2(-1)=2d0
      elseif ((case .eq. 'vlchwt') .or. (case .eq. 'vlchwn')
     .   .or. (case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
        taumin=0.0001d0
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0d0
        fx2(j)=0d0
        enddo
        fx1(2)=1d0
        fx2(-1)=1d0
      else
        write(6,*) 'Unimplemented process in lowint : case=',case
        stop 
      endif
      
      
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo 

      currentPDF=0

c--- do not calculate the flux if we're only checking the volume      
      if (case(1:4) .ne. 'vlch') then      
        flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      endif
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0d0
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
      
c--- calculate PDF's  
      if (((case .eq. 'qg_tbq') .or. (case .eq. '4ftwdk'))
     &     .and. (dynamicscale .eqv. .false.)) then
c--- for single top + b, make sure to use two different scales
        call fdist(ih1,xx(1),facscale_H,fx1_H)
        call fdist(ih2,xx(2),facscale_H,fx2_H)
        call fdist(ih1,xx(1),facscale_L,fx1_L)
        call fdist(ih2,xx(2),facscale_L,fx2_L)
      do j=-nf,nf
        if (j .eq. 0) then   ! heavy quark line has gluon init. state
          fx1(j)=fx1_H(j)
          fx2(j)=fx2_H(j)
        else
          fx1(j)=fx1_L(j)
          fx2(j)=fx2_L(j)
        endif
      enddo
      elseif (case(1:4) .ne. 'vlch') then
        if ((case .eq. 'bq_tpq') .or. (case .eq. 'qg_tbq')) then   
c--- single top: allow for different scales on each leg  
c---  (applies only if dynstring = 'DDIS')
          if (dynstring .eq. 'DDIS') then
            call fdist(ih1,xx(1),b1scale,fxb1)
            call fdist(ih2,xx(2),q2scale,fx2)
            call fdist(ih1,xx(1),q1scale,fx1)
            call fdist(ih2,xx(2),b2scale,fxb2)
        else          
            call fdist(ih1,xx(1),facscale,fx1)
            call fdist(ih2,xx(2),facscale,fx2)
        endif
        else   
c--- usual case            
          call fdist(ih1,xx(1),facscale,fx1)
          call fdist(ih2,xx(2),facscale,fx2)
        endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav    

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (omitgg) then 
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif

      if     ((case .eq. 'bq_tpq') .and. (dynstring .eq. 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (abs(j) .eq. 5) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (abs(k) .eq. 5) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0d0
      endif
      elseif ((case .eq. 'qg_tbq') .and. (dynstring .eq. 'DDIS')) then
c--- special case for dynamic scale in t-channel single top
        if     (j .eq. 0) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (k .eq. 0) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0d0
      endif
      else
c--- DEFAULT
        xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      endif

      xmsq=xmsq+xmsqjk
      xmsq_array(j,k)=xmsqjk
      
      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      if (currentPDF .eq. 0) then
        xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk
      endif
      
 20   continue
      enddo
      enddo

      if (currentPDF .eq. 0) then
        lowint=flux*pswt*xmsq/BrnRat
      endif
            
c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        PDFwgt(currentPDF)=flux*pswt*xmsq/BrnRat*wgt/itmx
        PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .     +PDFwgt(currentPDF)
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
      endif    

        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     .        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     .        )*wgt*flux*pswt/BrnRat/dfloat(itmx)

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      val=wgt*flux*pswt/BrnRat
      do j=-1,1
         do k=-1,1
!$omp atomic            
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       val*xmsq_bypart(j,k)
      enddo
      enddo

      val=lowint*wgt
      val2=val**2
      if (ieee_is_nan(val)) then
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
      if ((.not.unweight) .and. (dabs(val) .gt. wtmax)) then
!$omp critical(MaxWgt)
        wtmax=dabs(val)
!$omp end critical(MaxWgt)
      endif


      if (bin) then
         call nplotter(pjet,val,val2,0)
c---  POWHEG-style output if requested
         if (writepwg) then
!$omp critical(pwhgplotter)
            call pwhgplotter(p,pjet,val,0)
!$omp end critical(pwhgplotter)
         endif
      endif

c --- Check weights :
      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
        wtabs = dabs(val)
c--- note that r(ndim+2) is reserved for new_pspace in real, so unused at LO
        if (r(ndim+2) .lt. (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
          if (wtabs.lt.wtmax) then
            newwt = 1d0
          else
            newwt = wtabs/wtmax
          endif
          if (newwt .gt. 1.0d0) then
            write(6,*) 'WARNING : lowint : event with |weight| > 1.',
     &                 ' |weight| = ',newwt
          endif
c ---     just in case the weight was negative :
          newwt = newwt*dsign(1d0,val)
c         call nplotter(pjet,newwt,newwt,0)
!$omp critical(LowintWriteLHE)
          call mcfm_writelhe(pjet,xmsq_array,newwt)
!$omp end critical(LowintWriteLHE)
        endif
      endif

      return

 999  continue
      lowint=0d0
      ntotzero=ntotzero+1
      
      return
      end


