      function fragint(r,wgt)
      implicit none
      include 'types.f'
      real(dp):: fragint
           
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'debug.f'
      include 'vegas_common.f'
      include 'frag.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'efficiency.f'
      include 'noglue.f'
      include 'masses.f'
      include 'maxwt.f'
      include 'PDFerrors.f'
      include 'wts_bypart.f'
      include 'nflav.f'
      include 'ipsgen.f'
      include 'dm_params.f' 
      include 'outputflags.f' 
      include 'lastphot.f' 
      include 'x1x2.f'
      include 'bypart.f'
      include 'energy.f'
      include 'initialscales.f'
      integer:: ih1,ih2,j,k,sgnj,sgnk,nvec,pflav,pbarflav
      real(dp):: r(mxdim),wgt,pswt,
     & p(mxpart,4),flux,xmsq_bypart(-1:1,-1:1),
     & BrnRat,pjet(mxpart,4),val,val2,
     & xmsq,xmsqjk,W,msq(-nf:nf,-nf:nf),fx1(-nf:nf),fx2(-nf:nf),
     & ran2,msqdips(-nf:nf,-nf:nf),p_phys(mxpart,4),
     & wt34,wt345,wtprop,s34,s345,dot,wtips(4)
      real(dp):: m3,m4,m5
      logical:: bin,includedipole,vetow_2gam
      external qqb_w_g,qqb_z1jet,qqb_dirgam,qqb_2j_t,qqb_2j_s,
     & qqb_z2jetx,qqb_zaj,qqb_dm_monojet,qqb_gmgmjt,qqb_dirgam_g,
     & qqb_trigam_g
      common/density/ih1,ih2
      common/bin/bin
      common/BrnRat/BrnRat
      
c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2

!$omp atomic
      ntotshot=ntotshot+1
      fragint=0._dp

      W=sqrts**2

c----------------------------- GENERATE PHASE SPACE ---------------------------


c--- need to do something special for W_2gam and Z_2gam due to ipsgen
      if ((kcase==kW_2gam) .or. (kcase==kZ_2gam)) then
        npart=4
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,0,p,pswt,*999) 
        elseif (((ipsgen == 2) .and. (kcase==kZ_2gam)) .or.
     &          ((ipsgen == 3) .and. (kcase==kW_2gam))) then
            call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999) 
        else
           write(6,*) 'Parameter ipsgen not allowed'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif
        if (kcase==kW_2gam) then
c          if (vetow_2gam(p)) goto 999 ! partition PS according to ipsgen
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(3)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(3))
        endif
      else
c--- otherwise, use same PS generation as at LO
        call gen_lops(r,p,pswt,*999)
      endif
      z_frag=r(ndim)
      frag=.true.
                     
c--------------------------------- PHASE SPACE CUTS ---------------------------

      nvec=npart+2
      call dotem(nvec,p,s)
      
c----reject event if any s(i,j) is too small
c      call smalls(s,npart,*999)
c----reject event if any tau is too small
      call smalltau(p,npart,*999)

c--- the generated phase space point is the one that should be used
c--- in the calculation of the matrix elements;
c--- the physical momenta (p_phys) correspond to rescaling one of the photons
c--- "lastphot" by z_frag; this array should be used for cuts, plotting, etc.
      p_phys(:,:)=p(:,:)
      p_phys(lastphot,:)=z_frag*p(lastphot,:)

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p_phys) .eqv. .false.) then
        goto 999
      endif

c--- cut on z_frag      
      if((z_frag < 0.0001_dp) .or. (z_frag > 1._dp)) goto 999
       
      if (dynamicscale) then 
         call scaleset(initscale,initfacscale,p_phys)
      endif

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

c-------------------------- CALCULATE MATRIX ELEMENTS ---------------------------
      

c--------------------------------------------------------
c------------ NEED TO PASS P, P_PHYS TO ALL FRAGDIPS ROUTINES NOW
c---------------------------------------------------------
c-----------------------------------------------------------

      if     (kcase==kWgamma) then
         call qqb_wgam_frag(p,msq)
         call qqb_wgam_fragdips(p,p_phys,qqb_w_g,msqdips)
      elseif (kcase==kZgamma) then 
         call qqb_zgam_frag(p,msq)
         call qqb_zgam_fragdips(p,p_phys,qqb_z1jet,msqdips)
      elseif (kcase==kdirgam) then 
!         call qqb_dirgam_frag(p,msq)
!         call qqb_dirgam_fragdips(p,qqb_2j_t,qqb_2j_s,msqdips)
         msqdips(:,:)=0._dp 
!======== new format
         call qqb_dirgam_frag_combo(p,p_phys,msq) 
      elseif (kcase==kgamgam) then 
        call qqb_gamgam_frag(p,msq)
        call qqb_gamgam_fragdips(p,p_phys,qqb_dirgam,msqdips)
      elseif (kcase==ktrigam) then 
        call qqb_trigam_frag(p,msq)
        call qqb_trigam_fragdips(p,p_phys,qqb_gmgmjt,msqdips)
      elseif (kcase==kfourga) then 
         call qqb_fourgam_frag(p,msq) 
         call qqb_fourgam_fragdips(p,p_phys,qqb_trigam_g,msqdips)
      elseif (kcase==kgmgmjt) then 
         msqdips(:,:)=0._dp 
!====== new format 
        call qqb_gmgmjt_frag_combo(p,p_phys,msq)
!        call qqb_gmgmjt_fragdips(p,p_phys,msq,qqb_dirgam_g)
      elseif(kcase==kZ_2gam) then 
         call qqb_zaa_frag(p,msq) 
         call qqb_zaa_fragdips(p,p_phys,qqb_zaj,msqdips) 
      elseif(kcase==kZgajet) then 
         call qqb_zaj_frag(p,msq)
         call qqb_zaj_fragdips(p,p_phys,qqb_z2jetx,msqdips)
      elseif(kcase==kdm_gam) then 
         call qqb_dm_monophot_frag(p,msq) 
         call qqb_dm_monophot_fragdips(p,p_phys,qqb_dm_monojet,msqdips) 
      elseif(kcase==kW_2gam) then 
         stop
c         call qqb_waa_frag_combo(p,p_phys,msq) 
         msqdips(:,:)=0._dp 
      else
        write(6,*) 'Fragmentation MEs not available for this process.'
        write(6,*) 'kcase = ',kcase
        stop
      endif
      
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0._dp
      enddo
      enddo 

c--------------------------------------- INCLUDE PDF ---------------------------

    

      currentPDF=0

c--- do not calculate the flux if we're only checking the volume
      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0._dp
!      if (PDFerrors) then
!        call InitPDF(currentPDF)
!      endif
      
c--- calculate PDF's  
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

c--- sum of fragmentation contribution and integrated fragmentation dipoles
      xmsqjk=fx1(j)*fx2(k)*(msq(j,k)+msqdips(j,k))

      xmsq=xmsq+xmsqjk
      
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
        fragint=flux*pswt*xmsq/BrnRat
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

c      if (creatent) then
        wt_gg=xmsq_bypart(0,0)*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_gq=(xmsq_bypart(+1,0)+xmsq_bypart(-1,0)
     &        +xmsq_bypart(0,+1)+xmsq_bypart(0,-1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_qq=(xmsq_bypart(+1,+1)+xmsq_bypart(-1,-1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)
        wt_qqb=(xmsq_bypart(+1,-1)+xmsq_bypart(-1,+1)
     &        )*wgt*flux*pswt/BrnRat/real(itmx,dp)
c      endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      do j=-1,1
      do k=-1,1
!$omp atomic
        lord_bypart(j,k)=lord_bypart(j,k)+
     &       wgt*flux*pswt*xmsq_bypart(j,k)/BrnRat
      enddo
      enddo

      val=fragint*wgt
      val2=val**2
c--- update the maximum weight so far, if necessary
c---  but not if we are already unweighting ...
      if ((.not.unweight) .and. (abs(val) > wtmax)) then
!$omp critical(MaxWgt)
        wtmax=abs(val)
!$omp end critical(MaxWgt)
      endif

   

      if (bin) then
c ---   DSW. If the user has not selected to generate
c ---   events, still call nplotter here in order to
c ---   fill histograms/ntuples with weighted events :
        if (.not.evtgen) then
c!$omp critical(Plotter)
          call nplotter(pjet,val,val2,0)
c!$omp end critical(Plotter)
        endif
      endif

c --- Check weights :
c      if (unweight) then
c       write(6,*) 'Test weight ',val,' against max ',wtmax
c        wtabs = abs(val)
c        if (ran2() < (wtabs/wtmax)) then
c         write(6,*) 'Keep event with weight',val
c          if (wtabs<wtmax) then
c            newwt = 1._dp
c          else
c            newwt = wtabs/wtmax
c          endif
c          if (newwt > 1.0_dp) then
c            write(6,*) 'WARNING : fragint : event with |weight| > 1.',
c     +            ' |weight| = ',newwt
c          endif
c ---     just in case the weight was negative :
c          newwt = newwt*sign(1._dp,val)
!          call nplotter(pjet,newwt,newwt,0)
c ---     DSW. If I'm storing the event, I need to make a decision
c ---     about the flavours :
c          call decide_flavour(pflav,pbarflav)
c          call storeevent(pjet,newwt,pflav,pbarflav)
c        endif
c      endif

      return

 999  continue
      fragint=0._dp

      return
      end
      
