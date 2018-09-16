      !----------------------------------------------------------------------
      !      This replaces the original includedipole
      !      It calls the original and then the user one
      logical function includedipole(nd,ptrans)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      real(dp) ptrans(mxpart,4)
      integer nd
      logical mcfm_includedipole, userincludedipole
      
      ! it looks like all (nd .ne. 0) automatically have their momenta
      ! stored in the ptilde common; do this also for nd=0
      if (nd .eq. 0) call storeptilde(nd,ptrans)

      ! first call the original MCFM includedipole
      includedipole = mcfm_includedipole(nd,ptrans)

      ! then call the user one so that they can overturn the decision
      ! if they want
      includedipole = userincludedipole(nd,ptrans,includedipole)

      end 

      logical function mcfm_includedipole(nd,ptrans) 
     &                 result(mcfmincdipole)
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'notag.f'
      include 'taucut.f'
      include 'hdecaymode.f'
c---- SSbegin
      include 'reweight.f'
c---- SSend  

      real(dp) ptrans(mxpart,4),pjet(mxpart,4),rcut,pt,pttwo
      integer j,nd,isub
      logical gencuts,failedgencuts,photoncuts,makecuts,filterWbbmas,
     &     photonfailed,filterW_bjet,is_photon
      logical gencuts_VHbb
      integer count_photo,nphotons
      logical passed_frix,iso, passed_taucut
c      integer ij
c      real(dp) y32,y43,z3,z4,z5,z6
c      real(dp) dphizj,pt5sq,pt6sq,pt7sq

c      character*30 runstring
c      common/runstring/runstring
      common/rcut/rcut
      common/makecuts/makecuts
c---- SSbegin
c---- set default reweight to 1 (hence no reweighting)
      reweight = 1.0_dp
c---- SSend  

c--- default: include this contribution
      mcfmincdipole=.true.
      
c--- isub=1 for dipole subtractions, isub=0 for real radiation
      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

      nphotons=count_photo() 
      if (nphotons .gt. 0) then 
c--- Photons: Frixione isolation cuts if no fragmentation included 
         if (frag .eqv. .false.) then
            do j=3,mxpart 
               if(is_photon(j)) then 
                  call frix(ptrans,passed_frix,j,isub)
                  if(passed_frix.eqv..false.) then 
                     mcfmincdipole=.false.
                     return 
                  endif
               endif
            enddo
            call genclustphotons(ptrans,rcut,pjet,isub)
         else 
c--- Photons: not Frixione, need fragmentation and isolation
!---- do not want to cluster partons inside of jet cone, allow 
!---- isolation to describe these regions, therefore use 
!---- genclustphotons here 
            call genclustphotons(ptrans,rcut,pjet,isub)
!            call genclust2(ptrans,rcut,pjet,isub)
c---  Isolate photon
            do j=3,mxpart
            if (is_photon(j)) then 
               if (iso(ptrans,j,isub,nd) .eqv. .false.)then
                  mcfmincdipole=.false.
                  return 
               endif
            endif
            enddo
         endif
c--- check the photon cuts 
         photonfailed=photoncuts(pjet)
         if (photonfailed) then
            mcfmincdipole=.false.
            return
         endif 
      else
c--- No photons: the usual case
         call genclust2(ptrans,rcut,pjet,isub)
      endif


c--- perform mass cuts
      call masscuts(pjet,*999)
          
c--- fill ptilde array as persistent storage for the jet momenta
      ! GPS written in compact f90 form - in case we wish to copy it elsewhere
      ! (e.g. earlier, so that jets are defined even when includedipole is false)
      ptildejet(nd,1:npart+2,1:4)=pjet(1:npart+2,1:4)
      ! do j=1,4
      !   do i=1,npart+2
      !     ptildejet(nd,i,j)=pjet(i,j)
      !   enddo
      ! enddo
     
c--- for the Wbb process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing 
c--- NOTE: only for process numbers > 400 (20 and 25 should be handled normally) 
      if ((kcase==kWbbmas) .and. (nproc .gt. 400)) then
        mcfmincdipole=filterWbbmas()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif
      

c--- for the Wb+X process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing   
      if (kcase==kW_bjet) then
        mcfmincdipole=filterW_bjet()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif

     
      if (usescet) then
!==== special version for hadronic decays of Vector bosons at LO
         if(((kcase.eq.kWHbbar).or.(kcase.eq.kZHbbar)
     &        .or.(kcase.eq.kWH1jet).or.(kcase.eq.kZH1jet))
     &        .and.(hdecaymode=='bqba')) then 
            call maketaucut_bb(ptrans,jets,isub,passed_taucut)
         else
c--- branch for QT-cut 
!            passed_taucut=.true.
!            if ((npart==2).and.(isub==0)) then
!              continue
!            else
!              if (pttwo(3,4,ptrans) < taucut) passed_taucut=.false.
!            endif
!=== default
c---  for SCET calculation do not check jets, make tau cut instead
            call maketaucut(ptrans,pjet,jets,isub,passed_taucut)
         endif
        mcfmincdipole=passed_taucut
c           write(6,*) 'includedipole: nd,mcfmincdipole',nd,mcfmincdipole
c           if (passed_taucut .eqv. .false.) write(6,*) 'tau failed: ',nd
        if (mcfmincdipole .eqv. .false.) return
      else
c--- for a normal calculation,      
c--- if the number of jets is not correct, then do not include dipole
        if ((clustering .and. (jets .ne. nqcdjets-notag)
     &         .and. (inclusive .eqv. .false.)) .or.
     &      (clustering .and. (jets .lt. nqcdjets-notag)
     &         .and. (inclusive .eqv. .true.))) then
            mcfmincdipole=.false.
            return
        endif
      endif


!==== special case for WHbb and ZHbb, b's have special handling w.r.t tau cuts
!==== above, so call special cutting routine 
!      if((kcase==kWHbbar).or.(kcase.eq.kZHbbar)
!     & .or.(kcase==kWH1jet).or.(kcase.eq.kZH1jet)) then
!         if(usescet) then             
!            if(makecuts) then
!            failedgencuts=gencuts(pjet,jets)
!            if(failedgencuts) then             
!               mcfmincdipole=.false.
!               return
!            endif
!         endif
!         goto 99
!      endif
!      endif
      
c--- check the lepton cuts, if necessary
      if (makecuts) then
        failedgencuts=gencuts(pjet,jets)
        if (failedgencuts) then
          mcfmincdipole=.false.
          return
        endif
      endif
      
c      call Wgamcuts(ptrans,failedgencuts)
c      if (failedgencuts .eqv. .false.) mcfmincdipole=.false.

   99 continue
      
      return

  999 continue
      mcfmincdipole=.false.
      return    
      
      end
            




c--- This block of code may be reinstated above if required
cc--- special procedure for comparison with C. Oleari's
cc---  e+e- --> QQbg calculation; run Durham algorithm
c      if (runstring(1:5) .eq. 'carlo') then
c        call durhamalg(ptrans,npart-isub,y32,y43,z3,z4,z5,z6)
cc---     set up momentum array just in case it's needed elsewhere
c        do j=1,4
c          do i=1,npart+2
c          ptildejet(nd,i,j)=ptrans(i,j)
c          enddo
c        enddo
c      jets=1
c        if (y43 .gt. 0._dp) jets=2
c      if     ((kcase==kqq_tbg) .or. (kcase==kepem3j)) then
cc---     only keep events that have 3 jets when ycut=rcut
c        if ((y32 .lt. rcut) .or. (y43 .gt. rcut)) then
c          mcfmincdipole=.false.
c          jets=-1
c        endif
c        return
c      elseif (kcase==kqqtbgg) then
cc---     only keep events that have 4 jets when ycut=rcut
c        if (y43 .lt. rcut) then
c          mcfmincdipole=.false.
c          jets=-1
c        endif
c        return
c      else
c        write(6,*) 'Unexpected case in mcfmincdipole.f!'
c        stop
c      endif
c      endif
c

c--- This block of code may be reinstated above if required
cc--- added extra check here, to allow for analysis of G. Hesketh et al.
cc--- that requires Z+2 jets with only one jet within cuts, to obtain
cc--- prediction for Delta_phi(Z,jet) at NLO
c      if (runstring(1:6) .eq. 'dphizj') then
c        if    (nqcdjets .eq. 1) then
c        ij=5
c      elseif (nqcdjets .eq. 2) then
c        pt5sq=pjet(5,1)**2+pjet(5,2)**2
c        pt6sq=pjet(6,1)**2+pjet(6,2)**2
c        if (pt5sq .gt. pt6sq) then
c            ij=5
c        else
c          ij=6
c        endif
c      elseif (nqcdjets .eq. 3) then
c        pt5sq=pjet(5,1)**2+pjet(5,2)**2
c        pt6sq=pjet(6,1)**2+pjet(6,2)**2
c        pt7sq=pjet(7,1)**2+pjet(7,2)**2
c        if     (pt5sq .gt. max(pt6sq,pt7sq)) then
c            ij=5
c        elseif (pt6sq .gt. max(pt5sq,pt7sq)) then
c          ij=6
c        else
c          ij=7
c        endif
c      endif
c        dphizj=atan2(pjet(3,1)+pjet(4,1),pjet(3,2)+pjet(4,2))
c     .        -atan2(pjet(ij,1),pjet(ij,2))
c        if (dphizj .gt. pi) dphizj=twopi-dphizj
c        if (dphizj .lt. -pi) dphizj=twopi+dphizj
c      if (abs(dphizj) .gt. pi-1.e-1_dp) mcfmincdipole=.false.
c      endif
c      
