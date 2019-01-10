      logical function gencuts_input(pjet,njets)
************************************************************************
*   Author: J.M. Campbell, 5th December 2001                           *
*                                                                      *
*   This routine imposes a generic set of cuts that can be applied     *
*   to all processes in process.DAT, using the parton momenta in pjet  *
*   which have already passed through the jet clustering algorithm     *
*                                                                      *
*   Only a basic set of variables is tested:                           *
*     pt(lepton) > leptpt, eta(lepton) < leptrap, missing Et > misspt  *
*                                                                      *
*   If leptpt2 or leptrap2 is not zero, then any leptons beyond the    *
*   leading-pt one must instead satisfy:                               *
*     pt(lepton) > leptpt2, eta(lepton) < leptrap2                     *
*                                                                      *
*   For processes where one would like to apply jet-like cuts, but     *
*   no clustering has been performed, additional cuts apply:           *
*     pt(jet) > jetpt, eta(jet) < jetrap, R(jet1,jet2) > Rcut          *
*                                                                      *
*   Finally, if further process-specific cuts are necessary,           *
*   an appropriate second routine may be called                        *
*                                                                      *
*   Return TRUE if this point FAILS the cuts                           *
*                                                                      *
************************************************************************
      implicit none
      include 'bbproc.f'
      include 'constants.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'process.f'
      include 'plabel.f'
      include 'masses.f'
      include 'jetlabel.f'
      include 'notag.f'
      logical first,passedlept,is_lepton,is_photon,is_neutrino,
     & is_hadronic,failed
      integer njets,j,k,countb,bindex(mxpart),jindex,kindex,ib1,ib2
      integer countlept,leptindex(mxpart),countnu,
     & countjet,jetindex(mxpart),pntr,nuindex(mxpart)
      double precision pjet(mxpart,4),etvec(4),pZj(4),mZj
      double precision pt,etarap,etmiss,evtmisset,R,Rcut,etaj,etak,
     & etalept,mll,jetpt,jetrap,mnul,etabuffer
      logical hwwjetcuts
      double precision ht,qeta,mlbnu,merecon,reconcorr
      double precision dphi_ll,m_ll,mtrans,scut1,scut2
      double precision phill,phillcut,etajet2cut,mllcut,mnulcut
      double precision mZjcut,ptcheck,aetacheck
      integer ilep,igam,inu,countljet, countbjet,reconstr_top
      double precision pljet(mxpart,4),pbjet(mxpart,4)
      common/stopvars/ht,qeta,mlbnu,merecon,reconcorr
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
      common/rcut/Rcut
************************************************************************
*     Set-up the jet-like cut parameters here                          *
      parameter (jetpt=15d0,jetrap=2d0)
************************************************************************
      parameter(phillcut=1.2d0,etajet2cut=2.5d0,mllcut=15d0)
      parameter(mZjcut=40d0,mnulcut=30d0)
      data first/.true./
      save first,countlept,countnu,countb,leptindex,nuindex,bindex

      gencuts_input=.false.

      hwwjetcuts=.false.

      if (first) then
      first=.false.
c--- initialize counters and arrays that will be used to perform cuts
      countlept=0
      countnu=0
c--- lepton pt and rapidity cuts
      do j=3,mxpart
        if (is_lepton(j)) then
          countlept=countlept+1
          leptindex(countlept)=j
        endif
        if (is_neutrino(j)) then
           countnu=countnu+1
           nuindex(countnu)=j
        endif
      enddo
c--- Look for particles that should be treated as jets,
c--- so far only b decays from Z-bosons and
c--- hadronic decay of the W in diboson processes
c--- If these are present, we will do additional cuts
      countb=0
      do j=3,mxpart
         if ((plabel(j) .eq. 'qb') .or. (plabel(j) .eq. 'ab')
     .  .or. (plabel(j) .eq. 'qq') .or. (plabel(j) .eq. 'qa')) then
           countb=countb+1
           bindex(countb)=j
         endif
      enddo
c--- write-out the cuts we are using
      write(6,*)
      write(6,*)  '****************** Generic cuts ********************'
      write(6,*)  '*                                                  *'
      write(6,99) '*        pt(lepton)      >   ',leptpt,
     .                ' GeV            *'
      write(6,99) '*      |eta(lepton)|     <   ',leptrap,
     .                '                *'
      write(6,99) '*       pt(missing)      >   ',misspt,
     .                ' GeV            *'
      if ((leptpt2 .ne. 0d0) .or. (leptrap2 .ne. 0d0)) then
      write(6,99) '*     pt(2nd+ lepton)    >   ',leptpt2,
     .                ' GeV            *'
      write(6,99) '*   |eta(2nd+ lepton)|   <   ',leptrap2,
     .                '                *'
      endif
      if     (case.eq.'Wgamma') then
      if (mtrans34cut .lt. 0d0) then
      write(6,99) '*   (3,4,5) trans. mass  >   ',abs(mtrans34cut),
     .                ' GeV            *'
      else
      write(6,99) '* (e-gam,nu) trans. mass >   ',mtrans34cut,
     .                ' GeV            *'
      endif
      elseif (case.eq.'Zgamma') then
      write(6,99) '*    (3,4,5) inv. mass   >   ',mtrans34cut,
     .                ' GeV            *'
      elseif ( case .eq. 'Z_tjet') then
         write(6,99) '*      R(jet,lepton)     >   ',Rjlmin,
     .        '                *'
         write(6,99) '*    |mll-mz| <   ',mllcut,
     .                ' GeV            *'
         write(6,99) '*    |mZj-mt| <   ',mZjcut,
     .                ' GeV            *'

      else
      write(6,99) '*  (3,4) transverse mass >   ',mtrans34cut,
     .                ' GeV            *'
      endif
      write(6,99) '*      R(jet,lepton)     >   ',Rjlmin,
     .                '                *'
      write(6,99) '*     R(lepton,lepton)   >   ',Rllmin,
     .                '                *'
      write(6,99) '* |eta(jet1)-eta(jet2)|  >   ',delyjjmin,
     .                '                *'
      if (jetsopphem) then
      write(6,*) '*           eta(jet1) . eta(jet2)  <  0            *'
      endif
      if     (lbjscheme .eq. 1) then
      write(6,*) '*        eta(jet1)  <  eta(lept)  <  eta(jet2)     *'
      elseif (lbjscheme .ge. 2) then
      write(6,*) '*  eta(jet1)+Rcut  <  eta(lept)  <  eta(jet2)-Rcut *'
      lbjscheme=2
      endif
      if (countb .gt. 0) then
      write(6,*)  '*                                                  *'
      write(6,99) '*      pt(jet)       >   ',jetpt,
     .                ' GeV                *'
      write(6,99) '*    |eta(jet)|      <   ',jetrap,
     .                '                    *'
      write(6,99) '*   R(jet1,jet2)     >   ',Rcut,
     .                '                    *'
      endif
      if (hwwjetcuts) then
      write(6,99) '*   phi(lepton,lepton)   <   ',phillcut,
     .                '                *'
      write(6,99) '*  |eta(2nd jet)|        >   ',etajet2cut,
     .                '                *'
      write(6,99) '*  m(lepton,lepton)      <   ',mllcut,
     .                '                *'
      endif
      write(6,*)  '****************************************************'
      endif

C     Basic pt and rapidity cuts for lepton
      if     (countlept .eq. 1) then
          ptcheck=pt(leptindex(1),pjet)
          aetacheck=abs(etarap(leptindex(1),pjet))
          if ((ptcheck .lt. leptpt) .or. (aetacheck .gt. leptrap)) then
            gencuts_input=.true.
            return
          endif
          if ((aetacheck .gt. leptveto1min) .and.
     &        (aetacheck .lt. leptveto1max)) then
            gencuts_input=.true.
            return
          endif
      elseif (countlept .gt. 1) then
c--- loop over all the lepton possibilities for lepton 1 (j)
          j=0
  77      continue
          j=j+1
          passedlept=.true.
          ptcheck=pt(leptindex(j),pjet)
          aetacheck=abs(etarap(leptindex(j),pjet))
          if ((ptcheck .lt. leptpt) .or. (aetacheck .gt. leptrap)) then
            passedlept=.false.
            goto 78
          endif
          if ((aetacheck .gt. leptveto1min) .and.
     &        (aetacheck .lt. leptveto1max)) then
            passedlept=.false.
            goto 78
          endif
          do k=1,countlept
            if (k .ne. j) then
             ptcheck=pt(leptindex(k),pjet)
             aetacheck=abs(etarap(leptindex(k),pjet))
              if ((ptcheck.lt. leptpt2).or.(aetacheck.gt. leptrap2))then
                passedlept=.false.
              endif
              if ((aetacheck .gt. leptveto2min) .and.
     &            (aetacheck .lt. leptveto2max)) then
                passedlept=.false.
              endif
            endif
          enddo
  78      continue
c--- return to beginning if we failed and there are more leptons to try
          if ((passedlept .eqv. .false.).and.(j .lt. countlept)) goto 77
          gencuts_input=.not.(passedlept)
          if (gencuts_input) return
      endif

c--- missing energy cut
      evtmisset=etmiss(pjet,etvec)
      if ((evtmisset .lt. misspt) .and. (evtmisset .ne. 0d0)) then
        gencuts_input=.true.
        return
      endif

c--- mtrans34cut is used for three roles:
c---  1) Wgamma    --> mtrans34cut<0: transverse mass cut on (3,4,5) system
c---               --> mtrans34cut<0: transverse mass cut on (e-gam,nu) system
c---  2) Zgamma    --> invariant mass cut on (3,4,5) system
c---  3) otherwise --> transverse mass cut on (3,4) system
c---
      if (case.eq.'Wgamma') then
c--- cut on transverse mass of (3,4,5) system for Wgamma
        if (mtrans34cut .lt. 0d0) then
        mtrans=0d0
        do j=3,5
           mtrans=mtrans+dsqrt(pjet(j,1)**2+pjet(j,2)**2)
        enddo
        mtrans=mtrans**2
        do j=1,2
           mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2
        enddo
        mtrans=dsqrt(max(mtrans,0d0))
        if (mtrans .lt. abs(mtrans34cut)) then
           gencuts_input=.true.
           return
        endif
        else
c--- cut on (e-gam,nu) transverse mass for Wgamma
        if (is_neutrino(3)) then
          inu=3
          ilep=4
        else
          inu=4
          ilep=3
        endif
        igam=5
        mtrans=(pjet(ilep,4)+pjet(igam,4))**2
     &        -(pjet(ilep,1)+pjet(igam,1))**2
     &        -(pjet(ilep,2)+pjet(igam,2))**2
     &        -(pjet(ilep,3)+pjet(igam,3))**2
        mtrans=mtrans+(pjet(ilep,1)+pjet(igam,1))**2
     &               +(pjet(ilep,2)+pjet(igam,2))**2
        mtrans=dsqrt(max(mtrans,0d0))
     &        +dsqrt(pjet(inu,1)**2+pjet(inu,2)**2)
        mtrans=mtrans**2
        do j=1,2
           mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2
        enddo
        mtrans=dsqrt(max(mtrans,0d0))
        if (mtrans .lt. mtrans34cut) then
           gencuts_input=.true.
           return
        endif
        endif
c--- cut on invariant mass of (3,4,5) system for Zgamma
      elseif (case.eq.'Zgamma') then
        mtrans=(pjet(3,4)+pjet(4,4)+pjet(5,4))**2
        do j=1,3
          mtrans=mtrans-(pjet(3,j)+pjet(4,j)+pjet(5,j))**2
        enddo
        mtrans=dsqrt(max(mtrans,0d0))
        if(mtrans.lt.mtrans34cut) then
           gencuts_input=.true.
           return
        endif
      else
c--- cut on transverse mass of (3,4) pair otherwise
        mtrans=
     &   (pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2))
     &   /dsqrt((pjet(3,1)**2+pjet(3,2)**2)
     &         *(pjet(4,1)**2+pjet(4,2)**2))
        mtrans=2d0*dsqrt(pjet(3,1)**2+pjet(3,2)**2)
     &   *dsqrt(pjet(4,1)**2+pjet(4,2)**2)*(1d0-mtrans)
        mtrans=dsqrt(max(mtrans,0d0))
        if (mtrans .lt. mtrans34cut) then
          gencuts_input=.true.
        return
        endif
      endif

c--- lepton-lepton separation (if there are 2 or more leptons)
      if ((countlept .gt. 1)) then
        do j=1,countlept
        do k=j+1,countlept
          if (R(pjet,leptindex(j),leptindex(k)) .lt. Rllmin) then
            gencuts_input=.true.
            return
          endif
        enddo
        enddo
c--- extra cut on phi(lept,lept) for H(->WW)+jet search
        if (hwwjetcuts) then
          phill=
     .       (pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       +pjet(leptindex(1),2)*pjet(leptindex(2),2))
     .       /dsqrt((pjet(leptindex(1),1)**2+pjet(leptindex(1),2)**2)
     .             *(pjet(leptindex(2),1)**2+pjet(leptindex(2),2)**2))
          if (phill .lt. -0.999999999D0) phill=-1d0
          phill=dacos(phill)
          if (phill .gt. phillcut) then
            gencuts_input=.true.
            return
          endif
        endif
c--- extra cut on m(lept,lept) for H(->WW)+jet search
        if (hwwjetcuts) then
          mll=dsqrt(2d0*(
     .       +pjet(leptindex(1),4)*pjet(leptindex(2),4)
     .       -pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       -pjet(leptindex(1),2)*pjet(leptindex(2),2)
     .       -pjet(leptindex(1),3)*pjet(leptindex(2),3)))
          if (mll .gt. mllcut) then
            gencuts_input=.true.
            return
          endif
        endif
      endif

c -- mll lepton cut for Ztjet production
      if (case .eq. 'Z_tjet' .or. case .eq. 'Z_tdkj') then
         mll=dsqrt(2d0*(
     .        +pjet(leptindex(1),4)*pjet(leptindex(2),4)
     .        -pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .        -pjet(leptindex(1),2)*pjet(leptindex(2),2)
     .        -pjet(leptindex(1),3)*pjet(leptindex(2),3)))
c         if ( abs(mll-zmass) .ge. mllcut) then
c            gencuts_input=.true.
c            return
c         endif
      endif
      if (case .eq. 'Z_tdkj') then
         mnul=dsqrt(2d0*(
     .        +pjet(nuindex(1),4)*pjet(leptindex(3),4)
     .        -pjet(nuindex(1),1)*pjet(leptindex(3),1)
     .        -pjet(nuindex(1),2)*pjet(leptindex(3),2)
     .        -pjet(nuindex(1),3)*pjet(leptindex(3),3)))
c         if ( abs(mnul-wmass) .ge. mnulcut) then
c            gencuts_input=.true.
c            return
c         endif
      endif

c--- if there are no cuts on the jets - or no jets - we are done
      if ((Rjlmin .le. 0d0) .and. (delyjjmin .le. 0d0)) return
      if ((njets .eq. 0) .and. (countb .eq. 0)) return

c--- identify the jets
      countjet=0
      do j=3,mxpart
        if (is_hadronic(j)) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

c--- countjet will pick up the extra 'pp' needed for the real piece,
c--- therefore we should subtract 1 from this number
      if (countjet .gt. njets) countjet=countjet-1

      if ((njets .ne. countjet) .and. (notag .eq. 0)) then
        write(6,*) 'Something is wrong in gencuts_input.f -'
        write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
        stop
      endif

c -- identify b-jets
      call idjet(pjet,jetindex,countljet,countbjet,
     &     pljet,pbjet)

c--- extra cut on eta(2nd jet) for H(->WW)+jet search
      if ((hwwjetcuts) .and. (countjet .ge. 2)) then
        if (pt(jetindex(1),pjet) .gt. pt(jetindex(2),pjet)) then
          if (abs(etarap(jetindex(2),pjet)) .lt. etajet2cut) then
            gencuts_input=.true.
            return
          endif
        else
          if (abs(etarap(jetindex(1),pjet)) .lt. etajet2cut) then
            gencuts_input=.true.
            return
          endif
        endif
      endif
c--- jet-lepton separation (if there are 1 or more jets and leptons)
c -- applies to all jets
      if ((njets .gt. 0) .and. (countlept .gt. 0)) then
        do j=1,countlept
        do k=1,njets
           if (R(pjet,leptindex(j),jetindex(k)) .lt. Rjlmin) then
            gencuts_input=.true.
            return
          endif
        enddo
        enddo
      endif

c -- cut on Zj mass -- use light jet only
c -- only one light jet needs to be able to construct the top mass
      if (case .eq. 'Z_tjet' .or. case .eq. 'Z_tdkj') then
         reconstr_top=0
         do j=1,countljet
            pZj(:)=pljet(j,:)+pjet(leptindex(1),:)
     &           +pjet(leptindex(2),:)
            mZj=dsqrt(pZj(4)**2-pZj(1)**2-pZj(2)**2-pZj(3)**2)

            if ( abs(mZj-mt) .le. mZjcut) then
               reconstr_top=reconstr_top+1
            endif
         enddo
         if (reconstr_top .eq. 0) then
            gencuts_input=.true.
            return
         endif
      endif


c--- WBF-style cuts (if there are 2 or more jets)
      if ((njets .gt. 1)) then
c--- jet-jet rapidity separation
c--- j and k point to the two highest pt ('tagging') jets
        j=1
        k=2
        if (njets .eq. 3) then
          if     ( pt(jetindex(1),pjet) .lt.
     .      min(pt(jetindex(2),pjet),pt(jetindex(3),pjet)) ) then
            j=2
            k=3
          elseif ( pt(jetindex(2),pjet) .lt.
     .      min(pt(jetindex(1),pjet),pt(jetindex(3),pjet)) ) then
            j=1
            k=3
          endif
        endif
        if (abs(etarap(jetindex(j),pjet)-etarap(jetindex(k),pjet))
     .         .lt. delyjjmin) then
          gencuts_input=.true.
          return
        endif

c--- Requirement that jets be in opposite hemispheres
        if (jetsopphem) then
          if(etarap(jetindex(j),pjet)*etarap(jetindex(k),pjet) .gt. 0d0)
     .       then
            gencuts_input=.true.
            return
          endif
        endif

        if (lbjscheme .gt. 0) then
c--- Cut to require lepton to be between jets
          etabuffer=dble(lbjscheme-1)*Rcut
          etaj=etarap(jetindex(j),pjet)
          etak=etarap(jetindex(k),pjet)
          do pntr=1,countlept
            etalept=etarap(leptindex(pntr),pjet)
            if ( (etalept .lt. min(etaj,etak)+etabuffer) .or.
     .           (etalept .gt. max(etaj,etak)-etabuffer) ) then
              gencuts_input=.true.
              return
            endif
          enddo
        endif

      endif

c-- cuts on b-quarks
      if (bbproc) then
        call getbs(pjet,ib1,ib2)
        if ( (abs(etarap(ib1,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib1,pjet) .lt. ptbjetmin) ) gencuts_input=.true.
        if ( (abs(etarap(ib2,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib2,pjet) .lt. ptbjetmin) ) gencuts_input=.true.
      endif

c--- completed basic cuts
C--- if there are jet-like particles (see above), do more cuts
      if (countb .gt. 0) then
        do j=1,countb
          jindex=bindex(j)
          if (          (pt(jindex,pjet) .lt. jetpt) .or.
     .           (abs(etarap(jindex,pjet)) .gt. jetrap)) then
            gencuts_input=.true.
            return
          endif
          do k=j+1,countb
            kindex=bindex(k)
            if ((r(pjet,jindex,kindex) .lt. rcut)) then
              gencuts_input=.true.
              return
            endif
          enddo
        enddo
      endif

      return

   99 format(1x,a29,f6.2,a17)


      end




