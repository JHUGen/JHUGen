      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'phasemin.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'realwt.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'masses.f'
      include 'wts_bypart.f'
      include 'dipolescale.f'
      include 'stopscales.f'
      include 'flags.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'decay1q2a.f'
      include 'outputoptions.f'
      include 'breit.f'
      include 'dm_params.f' 
      include 'outputflags.f'
      include 'runstring.f'
      include 'bypart.f'
      include 'energy.f'
      include 'incldip.f'
      include 'nproc.f'
      include 'initialscales.f'

c--- APPLgrid - enable grids
c      include 'APPLinclude.f'
c      include 'qcdcouple.f'
c      double precision psCR
c--- APPLgrid - end

c---- SSbegin
      include 'reweight.f'
c---- SSend
cz
cz //
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec,ii,t
      double precision vector(mxdim),W,val,val2,valsum,xint,ptmp
      double precision fx1(-nf:nf),fx2(-nf:nf),
     . dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf),
     . fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqa(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4)
      double precision m3,m4,m5,R,Rbbmin
      double precision xmsq_bypart(0:maxd,-1:1,-1:1),xmsqjk,
     . plo(mxpart,4),pswtdip
      integer sgnj,sgnk
      common/xreal/xreal,xreal2
      common/Rbbmin/Rbbmin
      logical bin,failed
      logical includedipole,includereal
      double precision QandGint
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     . qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     . qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,
     . qqb_wh_ww,qqb_wh_ww_gs,
     . qqb_wh_zz,qqb_wh_zz_gs,
     . qqb_wbb_g,qqb_wbb_gs,
     . qqb_dirgam_g,qqb_dirgam_gs,qqb_hflgam_g,qqb_hflgam_gs,
     . qqb_trigam_g,qqb_trigam_gs,qqb_gmgmjt_g,qqb_gmgmjt_gs,
     . qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     . qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     . qqb_QQb_g,qqb_QQb_gs,
     . VV_Hqq_g,VV_Hqq_gs,VV_HWW_g,VV_HWW_gs,
     . gg_Hg,gg_H_gs,
     . gg_HWWgg,gg_HWWg_gs,gg_HZZgg,gg_HZZg_gs,
     . gg_Hgg,gg_Hg_gs,
     . gQ_zQ_g,gQ_zQ_gs,qqb_tbb_g,qqb_tbb_gs,
     . qqb_w_tndk_g,qqb_w_tndk_gs,
     . qqb_w_twdk_g,qqb_w_twdk_gs,qqb_w_twdk_gdk,qqb_w_twdk_gsdk,
     . qqb_zbjet_g,qqb_zbjet_gs,qqb_w_cjet_g,qqb_w_cjet_gs,
     . qqb_wbfromc_g,qqb_wbfromc_gs,
     . gg_hggg,gg_hgg_gs,qq_tchan_htqg,qq_tchan_htq_gs,
     . qg_tbq_g,qg_tbq_gs,qq_tbg_g,qq_tbg_gs,epem3j_g,epem3j_gs,
     . qq_tchan_ztqg_dk,qq_tchan_ztq_dk_gs,
     . qq_tchan_ztqg,qq_tchan_ztq_gs,
     . qqb_QQbdk_g,qqb_QQbdk_gs,qqb_gamgam_g,qqb_gamgam_gs,
     . qqb_zaj_gs,qqb_zaj_g,
     . qqb_gmgmjt_g_mad,
     . qqb_zaa_g,qqb_zaa_gs,
     . qqb_waa_g,qqb_waa_gs,
     . qqb_waa_g_mad,qqb_tottth_g,qqb_tottth_g_mad
      common/density/ih1,ih2
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
cz Add b fraction
      double precision bwgt
      common/btagging/ bwgt
      double precision msqtmp(0:maxd),bwgttmp(0:maxd)
      double precision realeventp(mxpart,4)
      common/realeventp/realeventp
      
      data bwgt / 0d0 /  ! in common block
cz // Add b fraction   Note: only msqtmp(0), bwgttmp(0) are used in nplotter.f
!      data p/56*0d0/
!$omp threadprivate(/pext/)


      QandGflag=.false.
      p(:,:)=0d0

      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0d0

      W=sqrts**2
      
c      if (first) then
c         write(6,*)
c         write(6,*) 'nmin=',nmin,',nmax=',nmax
c         write(6,*)
c         first=.false.
c      endif
      
c-- note: new_pspace now signifies multi-channel integration
      if (new_pspace) then
        call gen_lops(vector,plo,pswt,*999)
!        call writeout(plo)
        npart=npart+1
        call multichan(vector(ndim-2),vector(ndim-1),vector(ndim),
     &                 vector(ndim+2),plo,p,pswtdip,*999)
!        call writeout(p) 
!        pause
        pswt=pswt*pswtdip
      else
        call gen_realps(vector,p,pswt,*999)
      endif
            
      nvec=npart+2
      call dotem(nvec,p,s)
      
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 .gt.  1d0) .or. (xx2 .gt.  1d0)
     &.or.(xx1 .lt. xmin) .or. (xx2 .lt. xmin)) then
         goto 999
      endif

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
     
c--- extra cut to divide WQj/ZQj regions
      if ( (nproc .eq. 312) .or. (nproc .eq. 317)
     . .or.(nproc .eq. 322) .or. (nproc .eq. 327)
     . .or.(nproc .eq. 342) .or. (nproc .eq. 352)) then
        if (R(p,5,6) .lt. Rbbmin) goto 999
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 
 
      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
          msqLH(j,k)=0d0   ! for stop+b process
          msqHL(j,k)=0d0   ! for stop+b process
        enddo
        enddo
      endif
      
      do j=1,mxpart
      do k=1,4
        realeventp(j,k)=p(j,k)
      enddo
      enddo
      
c--- test to see whether we need Gflag and Qflag together
      if ( ((case .eq. 'W_2jet') .or. (case .eq. 'Z_2jet'))
     &.and. (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0d0
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif
      
c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue   
      
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,p)
        dipscale(0)=facscale
      endif
      
c---- generate collinear points that satisfy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c--- Calculate the required matrix elements
      if     (case .eq. 'W_only') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)
        call qqb_w_gs(p,msqc)
      elseif (case .eq. 'W_1jet') then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_w2jet(p,msq)
        call qqb_w1jet_gs(p,msqc)  
      elseif (case .eq. 'Wgamma') then
!         if(includereal) call singcheck(qqb_wgam_g,qqb_wgam_gs,p) ! Checked 08/27/02
        if (includereal) call qqb_wgam_g(p,msq)
        call qqb_wgam_gs(p,msqc)  
      elseif (case .eq. 'Wbfrmc') then
        if (includereal) call qqb_wbfromc_g(p,msq)
        call qqb_wbfromc_gs(p,msqc)  
      elseif (case .eq. 'W_cjet') then
c        call singcheck(qqb_w_cjet_g,qqb_w_cjet_gs,p) ! Checked 15/05/07
        if (includereal) call qqb_w_cjet_g(p,msq)
        call qqb_w_cjet_gs(p,msqc)
      elseif (case .eq. 'Zgamma') then
!         if(includereal) call singcheck(qqb_zgam_g,qqb_zgam_gs,p) ! Checked 01/04/11
        if (includereal) call qqb_zgam_g(p,msq)
        call qqb_zgam_gs(p,msqc)
      elseif (case .eq. 'Z_2gam') then
c         if(includereal) call singcheck(qqb_zaa_g,qqb_zaa_gs,p)
        if (includereal) call qqb_zaa_g(p,msq)
        call qqb_zaa_gs(p,msqc)
      elseif (case .eq. 'W_2gam') then
c        if(includereal) call singcheck(qqb_waa_g,qqb_waa_gs,p)
c        call compare_madgraph(p,qqb_waa_g,qqb_waa_g_mad)
         stop
c        if (includereal) call qqb_waa_g(p,msq)      
c        call qqb_waa_gs(p,msqc)
      elseif (case .eq. 'Zgajet') then
!        if(includereal) call singcheck(qqb_zaj_g,qqb_zaj_gs,p)     ! Checked 10/21/10
        if (includereal) call qqb_zaj_g(p,msq)      
        call qqb_zaj_gs(p,msqc)
      elseif (case .eq. 'Wbbmas') then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (case .eq. 'Wttmas') then
c        call singcheck(qqb_wbbm_g,qqb_wbbm_gs,p)     ! Checked 10/21/10
      if (includereal) call qqb_wbbm_g(p,msq)      
        call qqb_wbbm_gs(p,msqc)      
      elseif (case .eq. 'Wbbbar') then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (case .eq. 'W_2jet') then      
c        call singcheck(qqb_w2jet_g,qqb_w2jet_gs,p) ! Re-checked June 09
        if (includereal)  call qqb_w2jet_g(p,msq)
        call qqb_w2jet_gs(p,msqc)
      elseif (case .eq. 'Z_only') then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'Z_1jet') then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        if (includereal) call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (case .eq. 'Z_2jet') then
c        call singcheck(qqb_z2jet_g,qqb_z2jet_gs,p) ! Re-checked June 09
        if (includereal) call qqb_z2jet_g(p,msq)  
        call qqb_z2jet_gs(p,msqc) 
      elseif (case .eq. 'Zbbbar') then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        if (includereal) call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WWqqbr') then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        if (includereal) call qqb_ww_g(p,msq)      
        call qqb_ww_gs(p,msqc)      
      elseif (case .eq. 'WWqqdk') then
        if (includereal) call dkqqb_ww_g(p,msq)      
        call dkqqb_ww_gs(p,msqc)      
      elseif (case .eq. 'WZbbar') then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_wz_g(p,msq)      
        call qqb_wz_gs(p,msqc)      
      elseif (case .eq. 'ZZlept') then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        if (includereal) call qqb_zz_g(p,msq)      
        call qqb_zz_gs(p,msqc)      
      elseif (case .eq. 'WHbbar') then
c        call singcheck(qqb_wh_g,qqb_wh_gs,p)
        if (includereal) call qqb_wh_g(p,msq)      
        call qqb_wh_gs(p,msqc)     
      elseif (case .eq. 'WH__WW') then
c        call singcheck(qqb_wh_ww_g,qqb_wh_ww_gs,p)
        if (includereal) call qqb_wh_ww_g(p,msq)      
        call qqb_wh_ww_gs(p,msqc)
      elseif (case .eq. 'WH__ZZ') then
c        call singcheck(qqb_wh_zz_g,qqb_wh_zz_gs,p)
        if (includereal) call qqb_wh_zz_g(p,msq)      
        call qqb_wh_zz_gs(p,msqc)   
      elseif (case .eq. 'WHgaga') then
c        call singcheck(qqb_wh_gaga_g,qqb_wh_gaga_gs,p)
        if (includereal) call qqb_wh_gaga_g(p,msq)      
        call qqb_wh_gaga_gs(p,msqc) 
      elseif (case .eq. 'ZHbbar') then
c        call singcheck(qqb_zh_g,qqb_zh_gs,p)
        if (includereal) call qqb_zh_g(p,msq)      
        call qqb_zh_gs(p,msqc)     
      elseif (case .eq. 'ZHgaga') then
c        call singcheck(qqb_zh_gaga_g,qqb_zh_gaga_gs,p)
        if (includereal) call qqb_zh_gaga_g(p,msq)      
        call qqb_zh_gaga_gs(p,msqc)     
      elseif (case .eq. 'ZH__WW') then
c        call singcheck(qqb_zh_ww_g,qqb_zh_ww_gs,p)
        if (includereal) call qqb_zh_ww_g(p,msq)      
        call qqb_zh_ww_gs(p,msqc)     
      elseif (case .eq. 'ZH__ZZ') then
c        call singcheck(qqb_zh_zz_g,qqb_zh_zz_gs,p)
        if (includereal) call qqb_zh_zz_g(p,msq)      
        call qqb_zh_zz_gs(p,msqc)     
      elseif (case .eq. 'dirgam') then
!        if (includereal) call singcheck(qqb_dirgam_g,qqb_dirgam_gs,p) 
        if (includereal) call qqb_dirgam_g(p,msq)
        call qqb_dirgam_gs(p,msqc)
      elseif (case .eq. 'hflgam') then
c        if (includereal) call singcheck(qqb_hflgam_g,qqb_hflgam_gs,p)
        if (includereal) call qqb_hflgam_g(p,msq)
        call qqb_hflgam_gs(p,msqc)
      elseif (case .eq. 'gamgam') then
c        if (includereal) call singcheck(qqb_gamgam_g,qqb_gamgam_gs,p)
        if (includereal) call qqb_gamgam_g(p,msq)
        call qqb_gamgam_gs(p,msqc)
      elseif (case .eq. 'gmgmjt') then
c        if (includereal) call singcheck(qqb_gmgmjt_g,qqb_gmgmjt_gs,p)
        if (includereal) call qqb_gmgmjt_g(p,msq)
        call qqb_gmgmjt_gs(p,msqc)
      elseif (case .eq. 'trigam') then
c        if (includereal) call singcheck(qqb_trigam_g,qqb_trigam_gs,p)
        if (includereal) call qqb_trigam_g(p,msq)
        call qqb_trigam_gs(p,msqc)
      elseif (case.eq.'fourga') then 
c        if (includereal) call singcheck(qqb_fourgam_g,qqb_fourgam_gs,p)       
        if (includereal) call qqb_fourgam_g(p,msq)
        call qqb_fourgam_gs(p,msqc)
       elseif (case .eq. 'ggfus0') then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         if (includereal) call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
       elseif (case .eq. 'Higaga') then
c         call singcheck(gg_hgamgamg,gg_hgamgam_gs,p)
         if (includereal) call gg_hgamgamg(p,msq)
         call gg_hgamgam_gs(p,msqc)
       elseif (case .eq. 'Hi_Zga') then
c         call singcheck(gg_hzgamg,gg_hzgam_gs,p)
         if (includereal) call gg_hzgamg(p,msq)
         call gg_hzgam_gs(p,msqc)
      elseif ((case .eq. 'HWW_4l') .or. (case .eq. 'HWW2lq')) then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        if (includereal) call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)      
      elseif (case .eq. 'HWWdkW') then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (case .eq. 'HWWdkW') then
        if (includereal) call dkqqb_hww_g(p,msq)      
        call dkqqb_hww_gs(p,msqc)      
      elseif (case .eq. 'HZZ_4l') then
c        call singcheck(qqb_hzz_g,qqb_hzz_gs,p)
        if (includereal) call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (case .eq. 'H_1jet') then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        if (includereal) call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
c         call singcheck(qqb_QQbdk_g,qqb_QQbdk_gs,p) ! Checked 15/8/08
        if (includereal) call qqb_QQbdk_g(p,msq)  
        call qqb_QQbdk_gs(p,msqc) 
      elseif ((case .eq. 'tt_ldk') .or. (case .eq. 'tt_hdk')) then
        if     (decay1q2a .eq. 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_QQb_g(p,msq)
          call dk1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a .eq. 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_QQb_g(p,msq)
          call dk2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (case .eq. 'tt_udk') then
        if     (decay1q2a .eq. 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1uqqb_QQb_g(p,msq)
          call dk1uqqb_QQb_gs(p,msqc)
      elseif (decay1q2a .eq. 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2uqqb_QQb_g(p,msq)
          call dk2uqqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (case .eq. 'tthWdk') then
        if     (decay1q2a .eq. 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dkW1qqb_QQb_g(p,msq)
          call dkW1qqb_QQb_gs(p,msqc)
      elseif (decay1q2a .eq. 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dkW2qqb_QQb_g(p,msq)
          call dkW2qqb_QQb_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (case .eq. 'tt_bbu') then
c         call singcheck(qqb_QQbdku_g,qqb_QQbdku_gs,p) !
c         pause
        if (includereal) call qqb_QQbdku_g(p,msq)  
        call qqb_QQbdku_gs(p,msqc) 
      elseif ((case .eq. 'tt_tot') .or. (case .eq. 'cc_tot')
     .   .or. (case .eq. 'bb_tot')) then
c        call singcheck(qqb_QQb_g,qqb_QQb_gs,p)
        if (includereal) call qqb_QQb_g(p,msq)
        call qqb_QQb_gs(p,msqc)
      elseif (case .eq. 'bq_tpq') then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbb_g(p,msq)
        call qqb_tbb_gs(p,msqc)
      elseif (case .eq. 't_bbar') then
c        call singcheck(qqb_tbb_g,qqb_tbb_gs,p)
        if (includereal) call qqb_tbbdk_g(p,msq)
        call qqb_tbbdk_gs(p,msqc)
      elseif (case .eq. 'ttdkay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
        if (includereal) call bq_tpq_gdk(p,msq)
        call bq_tpq_gsdk(p,msqc)
      elseif (case .eq. 'tdecay') then
c        call singcheck(qqb_tbb_g_new,qqb_tbb_gs,p)
      if (includereal) call dkqqb_tbbdk_g(p,msq)
        call dkqqb_tbbdk_gs(p,msqc)       
       elseif (case .eq. 'W_tndk') then
c        call singcheck(qqb_w_tndk_g,qqb_w_tndk_gs,p)      ! Checked 12/3/04
        if (includereal) call qqb_w_tndk_g(p,msq)
        call qqb_w_tndk_gs(p,msqc)
      elseif (case .eq. 'W_twdk') then
c        call singcheck(qqb_w_twdk_g,qqb_w_twdk_gs,p)            ! Checked 2/4/05
        if (includereal) call qqb_w_twdk_g(p,msq)
        call qqb_w_twdk_gs(p,msqc)
      elseif (case .eq. 'Wtdkay') then
        if (includereal) call dkqqb_w_twdk_g(p,msq)
        call dkqqb_w_twdk_gs(p,msqc)
      elseif ( (case .eq. 'qq_ttw')) then 
        if (includereal) call qqb_ttw_g(p,msq)
        call qqb_ttw_gs(p,msqc)
      elseif ( (case .eq. 'ttwldk')) then
        if     (decay1q2a .eq. 1) then
c------ radiation in the decay of the top quark
        if (includereal) call dk1qqb_ttw_g(p,msq)
          call dk1qqb_ttw_gs(p,msqc)
      elseif (decay1q2a .eq. 2) then 
c------ radiation in the decay of the anti-top quark
        if (includereal) call dk2qqb_ttw_g(p,msq)
          call dk2qqb_ttw_gs(p,msqc)
      else
        write(6,*) 'Problem generating PS, decay1q2a=',decay1q2a
        stop
      endif
      elseif (case .eq. 'ggfus1') then
c        call singcheck(gg_hgg,gg_hg_gs,p)
        if (includereal) call gg_hgg(p,msq)
        call gg_hg_gs(p,msqc)
      elseif (case .eq. 'Hgagaj') then
c        call singcheck(gg_hgagagg,gg_hgagag_gs,p)
        if (includereal) call gg_hgagagg(p,msq)
        call gg_hgagag_gs(p,msqc)
      elseif (case .eq. 'HWWjet') then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWgg(p,msq)
        call gg_hWWg_gs(p,msqc)
      elseif (case .eq. 'HWW2jt') then
c        call singcheck(gg_hWWgg,gg_hWWg_gs,p)
        if (includereal) call gg_hWWggg(p,msq)
        call gg_hWWgg_gs(p,msqc)
      elseif (case .eq. 'HZZjet') then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZgg(p,msq)
        call gg_hZZg_gs(p,msqc)
      elseif (case .eq. 'HZZ2jt') then
c        call singcheck(gg_hZZgg,gg_hZZg_gs,p)
        if (includereal) call gg_hZZggg(p,msq)
        call gg_hZZgg_gs(p,msqc)
c      write(6,*) msq
c      pause
      elseif (case .eq. 'qq_Hqq') then
c        call singcheck(VV_Hqq_g,VV_Hqq_gs,p)
        if (includereal) call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
      elseif (case .eq. 'qq_Hgg') then
c        call singcheck(VV_Hgaga_g,VV_Hgaga_gs,p)
        if (includereal) call VV_Hgaga_g(p,msq)
        call VV_Hgaga_gs(p,msqc)
      elseif (case .eq. 'qq_HWW') then
c        call singcheck(VV_HWW_g,VV_HWW_gs,p)
        if (includereal) call VV_HWW_g(p,msq)
        call VV_HWW_gs(p,msqc)
      elseif (case .eq. 'qq_HZZ') then
c        call singcheck(VV_HZZ_g,VV_HZZ_gs,p)
        if (includereal) call VV_HZZ_g(p,msq)
        call VV_HZZ_gs(p,msqc)
      elseif (case .eq. 'ggfus2') then
c        call singcheck(gg_hggg,gg_hgg_gs,p)            ! Checked 10/29/09
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (case .eq. 'gagajj') then
c        call singcheck(gg_hggg,gg_hgg_gs,p)
        if (includereal) call gg_hggg(p,msq)
        call gg_hgg_gs(p,msqc)
      elseif (case .eq. 'qg_tbq') then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)       ! Checked 2/4/08
       if (includereal) call qg_tbq_g(p,msq)
       call qg_tbq_gs(p,msqc)
      elseif (case .eq. '4ftwdk') then
c        call singcheck(qg_tbq_g,qg_tbq_gs,p)
       if (includereal) call qg_tbqdk_g(p,msq)
       call qg_tbqdk_gs(p,msqc)
      elseif (case .eq. 'dk_4ft') then
c       if (includereal) call dkqg_tbqdk_g_old(p,msq)  
       if (includereal) call dkqg_tbqdk_g(p,msq)  
       call dkqg_tbqdk_gs(p,msqc)
      elseif (case .eq. 'qq_tbg') then
c        call singcheck(qq_tbg_g,qq_tbg_gs,p)       ! Checked 8/9/08
       if (includereal) call qq_tbg_g(p,msq)
       call qq_tbg_gs(p,msqc)
      elseif (case .eq. 'epem3j') then
c        call singcheck(epem3j_g,epem3j_gs,p)       ! Checked 17/11/08
       if (includereal) call epem3j_g(p,msq)
       call epem3j_gs(p,msqc)
      elseif (case .eq. 'gQ__ZQ') then
c        call singcheck(gQ_zQ_g,gQ_zQ_gs,p)
        if (includereal) call gQ_zQ_g(p,msq)
        call gQ_zQ_gs(p,msqc)
      elseif (case .eq. 'Z_bjet') then
c        call singcheck(qqb_zbjet_g,qqb_zbjet_gs,p)      ! Checked 07/18/05
        if (includereal) call qqb_zbjet_g(p,msq)
        call qqb_zbjet_gs(p,msqc)
      elseif (case .eq. 'W_bjet') then
c        call singcheck(qqb_wbjet_g,qqb_wbjet_gs,p) ! Rechecked 14/3/08
        if (includereal) call qqb_wbjet_g(p,msq)
        call qqb_wbjet_gs(p,msqc)
      elseif (case .eq. 'Wcsbar') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        if (includereal) call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'Wcs_ms') then
        if (includereal) call qqb_w_cjet(p,msq)
        ndmax=0
      elseif (case .eq. 'Z_tjet') then
c        if (includereal) call singcheck(qq_tchan_ztqg,qq_tchan_ztq_gs,p)
        if (includereal) call qq_tchan_ztqg(p,msq)
        call qq_tchan_ztq_gs(p,msqc)
      elseif (case .eq. 'Z_tdkj') then
c        if (includereal) call singcheck(qq_tchan_ztqg_dk,
c     &                                  qq_tchan_ztq_dk_gs,p)
        if (includereal) call qq_tchan_ztqg_dk(p,msq)
        call qq_tchan_ztq_dk_gs(p,msqc)
      elseif (case .eq. 'H_tjet') then
c        if (includereal) call singcheck(qq_tchan_htqg,qq_tchan_htq_gs,p)
        if (includereal) call qq_tchan_htqg(p,msq)
        call qq_tchan_htq_gs(p,msqc)
      elseif (case .eq. 'H_tdkj') then
        if (includereal) call qq_tchan_htqg_dk(p,msq)
        call qq_tchan_htq_dk_gs(p,msqc)
      elseif (case .eq. 'tottth') then
c        if (includereal) call qqb_tottth_g(p,msq)
c        call qqb_tottth_gs(p,msqsc)
c        call singcheck(qqb_tottth_g,qqb_tottth_gs,p)
c        call compare_madgraph(p,qqb_tottth_g,qqb_tottth_g_mad)
      elseif (case .eq. 'dm_jet') then 
         if(includereal) call qqb_dm_monojet_g(p,msq)       
         call qqb_dm_monojet_gs(p,msqc)
      elseif (case.eq.'dm_gam') then
!         if(includereal) call singcheck(qqb_dm_monophot_g
!     &        ,qqb_dm_monophot_gs,p)
         if(includereal) call qqb_dm_monophot_g(p,msq) 
         call qqb_dm_monophot_gs(p,msqc)

      endif
      
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      do j=-1,1
      do k=-1,1
      xmsq_bypart(nd,j,k)=0d0
      enddo
      enddo
      enddo
      
      currentPDF=0
            
      flux=fbGeV2/(two*xx1*xx2*W)
c--- for mlm study, divide by (Ecm)**2=W
c      if (runstring(1:3) .eq. 'mlm') then
c      flux=flux/W
c      endif

c--- initialize a PDF set here, if calculating errors
  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
cz
      msqtmp(nd)=0d0
      bwgttmp(nd)=0d0
cz //
      enddo
      if (PDFerrors) then
        call InitPDF(currentPDF)
      endif
         
c--- calculate PDF's  
      if (dynamicscale) then
        do nd=ndmax,0,-1  ! so that fx1,fx2 correct for real kinematics
          if (dipscale(nd) .lt. 1d-8) then        
c--- in case dipole is not used, set up dummy value of scale for safety
c--- and set all PDF entries to zero
          dipscale(nd)=dipscale(0)
          do j=-nf,nf
            fx1(j)=0d0
            fx2(j)=0d0
          enddo
        else
            call fdist(ih1,xx1,dipscale(nd),fx1)
            call fdist(ih2,xx2,dipscale(nd),fx2)
            do j=-nf,nf
            dipfx1(nd,j)=fx1(j)
            dipfx2(nd,j)=fx2(j)
          enddo
        endif
      enddo
        if ((case .eq. 'qg_tbq') .or. (case .eq. '4ftwdk')) then
        do j=-nf,nf
        fx1_H(j)=fx1(j)
        fx1_L(j)=fx1(j)
        fx2_H(j)=fx2(j)
        fx2_L(j)=fx2(j)
        enddo
      endif
      else
        if ((case .eq. 'qg_tbq') .or. (case .eq. '4ftwdk')) then
c--- for single top + b, make sure to use two different scales
          call fdist(ih1,xx1,facscale_H,fx1_H)
          call fdist(ih2,xx2,facscale_H,fx2_H)
          call fdist(ih1,xx1,facscale_L,fx1_L)
          call fdist(ih2,xx2,facscale_L,fx2_L)
        do j=-nf,nf
          if (j .eq. 0) then  ! heavy quark line has gluon init. state
            fx1(j)=fx1_H(j)
            fx2(j)=fx2_H(j)
          else
            fx1(j)=fx1_L(j)
            fx2(j)=fx2_L(j)
          endif
        enddo
        else
c--- for comparison with C. Oleari's e+e- --> QQbg calculation
c            if (runstring(1:5) .eq. 'carlo') then
c            flux=1d0/2d0/W/(as/twopi)**2
c--- divide out by (ason2pi) and then the "LO" massless DY process
c          flux=flux/(aveqq*xn*fourpi*(gwsq/fourpi)**2/3d0/sqrts**2)
c            flux=flux/(xn/8d0)
c          do j=-nf,nf
c          fx1(j)=0d0
c          fx2(j)=0d0
c          enddo
c          fx1(0)=1d0
c          fx1(1)=1d0
c          fx2(0)=1d0
c          fx2(1)=1d0
c          else   
c--- usual case            
            call fdist(ih1,xx1,facscale,fx1)
            call fdist(ih2,xx2,facscale,fx2)
c        endif
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

      if ((case .eq. 'Wcsbar').and.(j .ne. 4).and.(k .ne. 4)) goto 20

      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
         if (dynamicscale) then         
             xmsq(nd)=xmsq(nd)+dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
         else
             xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         endif
         enddo
      else

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

c--- for single top + b, make sure to use two different scales
         if ((case .eq. 'qg_tbq') .or. (case .eq. '4ftwdk')) then
           xmsqjk=fx1_L(j)*fx2_H(k)*msqLH(j,k)
     .           +fx1_H(j)*fx2_L(k)*msqHL(j,k)
         else
c--- usual case            
           xmsqjk=fx1(j)*fx2(k)*msq(j,k)
       endif
       
         xmsq(0)=xmsq(0)+xmsqjk
cz
cz Extract fraction with b in final state, store in common
cz
c     nproc=161 then t-channel: t
c     isub = 1, nwz = 1
         if (nproc.eq.161) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.lt.0).or.((k.gt.0).and.
     &           (k.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.lt.0).or.((j.gt.0).and.
     &           (j.ne.5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
c     nproc=166 then t-channel: t~
         if (nproc.eq.166) then ! t
            msqtmp(0)=msqtmp(0)+xmsqjk
            if ((j.eq.0).and.((k.gt.0).or.((k.lt.0).and.
     &           (k.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
            if ((k.eq.0).and.((j.gt.0).or.((j.lt.0).and.
     &           (j.ne.-5)))) then
               bwgttmp(0)=bwgttmp(0)+xmsqjk
            endif
         endif
cz // end fill index 0

         if (currentPDF .eq. 0) then
           xmsq_bypart(0,sgnj,sgnk)=xmsq_bypart(0,sgnj,sgnk)+xmsqjk
         endif
         do nd=1,ndmax
         if (dynamicscale) then         
             xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
         else
             xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
         endif
           xmsq(nd)=xmsq(nd)+xmsqjk
           if (currentPDF .eq. 0) then
             xmsq_bypart(nd,sgnj,sgnk)=xmsq_bypart(nd,sgnj,sgnk)+xmsqjk
           endif
         enddo
         
      endif
 20   continue

      enddo
      enddo

c--- loop over all PDF error sets, if necessary
      if (PDFerrors) then
        do nd=0,ndmax 
          PDFxsec_nd(currentPDF,nd)=xmsq(nd)
        enddo
        currentPDF=currentPDF+1
        if (currentPDF .le. maxPDFsets) goto 777
c--- reset xmsq to the central PDF values
        do nd=0,ndmax 
          xmsq(nd)=PDFxsec_nd(0,nd)
        enddo
      endif    

      realint=0d0
      xint=0d0
cz
      bwgt=0d0
cz //

      valsum=0d0 ! running total of weights at this point

c--- zero out temporary histograms
c      if (bin) call zerorealhistos
      if (bin) call smartzero

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
c        if (creatent) then
          wt_gg=xmsq_bypart(nd,0,0)*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_gq=(xmsq_bypart(nd,+1,0)+xmsq_bypart(nd,-1,0)
     .          +xmsq_bypart(nd,0,+1)+xmsq_bypart(nd,0,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qq=(xmsq_bypart(nd,+1,+1)+xmsq_bypart(nd,-1,-1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
          wt_qqb=(xmsq_bypart(nd,+1,-1)+xmsq_bypart(nd,-1,+1)
     .          )*wgt*flux*pswt/BrnRat/dfloat(itmx)
c        endif
        failed=.false.
        
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo

          if (incldip(nd)) incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997
        endif

c---if it does, add to total
!        xint=xint+xmsq(nd)
c---  SSbegin
        xint=xint+xmsq(nd)*reweight
c---  SSend

        val=wgt*flux*pswt/BrnRat
        do j=-1,1
        do k=-1,1
!$omp atomic
          lord_bypart(j,k)=lord_bypart(j,k)+
     .         val*xmsq_bypart(nd,j,k)
        enddo
        enddo

        val=xmsq(nd)*wgt
        val2=val**2
      
      valsum=valsum+val
      
cz Fill bwgt if needed
        if(dabs(msqtmp(nd)).gt.0d0) bwgt=bwgttmp(nd)/msqtmp(nd)
cz //
        
c--- update PDF errors
        if (PDFerrors) then
          do currentPDF=0,maxPDFsets        
          PDFwgt(currentPDF)=
     .       flux*pswt*PDFxsec_nd(currentPDF,nd)/BrnRat*wgt/itmx
          PDFxsec(currentPDF)=PDFxsec(currentPDF)
     .       +PDFwgt(currentPDF)
          enddo           
        endif
                
c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)
          call nplotter(pjet,val,val2,nd)
c--- POWHEG-style output if requested
          if (writepwg) then
            if (nd .eq. 0) then 
!$omp critical(pwhgplotter)
              call pwhgplotter(p,pjet,val,nd)
!$omp end critical(pwhgplotter)
            else
             do j=1,mxpart
              do k=1,4
              q(j,k)=ptilde(nd,j,k)
              enddo
              enddo
!$omp critical(pwhgplotter)
              call pwhgplotter(q,pjet,val,nd)
!$omp end critical(pwhgplotter)
            endif
          endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

c 998  continue

c--- add temporary histograms to cumulative totals
      if (bin) then
         call smartadd(wgt) 
      endif
c--- update the maximum weight so far, if necessary
      if (dabs(valsum) .gt. wtmax) then
!$omp critical(MaxWgt)
        wtmax=dabs(valsum)
!$omp end critical(MaxWgt)
      endif
      
      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      
c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+realint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag), calling includedipole first to reset
c--- the values of "jets" and "jetlabel"
          includereal=includedipole(0,p)
          Qflag=.true.
          Gflag=.false.
          goto 44
        else
c--- return both to .true. and assign value to realint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          realint=QandGint
        endif
      endif

      return

 999  realint=0d0
      ntotzero=ntotzero+1
c--- safety catch
      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif      
      return
      end

