      function gencuts_WZjj(pjet,njets)
       implicit none
      include 'types.f'
      logical:: gencuts_WZjj
c -- these cuts taken from CMS FCNC search, arxiv:1208.0957
c -- R. Rontsch 2103-03-01

      include 'leptcuts.f'
      include 'jetcuts.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'plabel.f'
      include 'masses.f'
      include 'nwz.f'
      include 'runstring.f'
      include 'notag.f'
      integer:: njets
      real(dp):: pjet(mxpart,4)
      integer:: leptindex(3), jetindex(mxpart),countjet
      integer:: j,k,jj,kk, Zpair1(2),Zpair2(2),Wpair1(2),Wpair2(2),
     &     Wpair(2),Zpair(2),nuindex
      logical:: first, diffflav
      real(dp):: mll1,mll2,threemass,twomass,
     & r,pt,etarap
      real(dp):: mllmax,mllmin,Rleptisol,Zisol,Wisol,
     &     mZjcut,mWbcut,incone(3),mWbmin,mWbmax,
     &     mZjmin,mZjmax,STmin,ST,isol(3),isol1(3),isol2(3),
     &     mZj1,mZj2,mWj1,mWj2
      parameter (mllmin=60._dp,mllmax=120._dp,Rleptisol=0.3_dp,Zisol=0.125_dp,
     &     Wisol=0.1_dp,mZjcut=25._dp,mWbcut=35._dp,mZjmin=100._dp,
     &     mZjmax=250._dp,mWbmin=100._dp,mWbmax=250._dp,STmin=250._dp)
      character*4 cut_id
      data first/.true./

      gencuts_WZjj=.false.
c -- allows for same flavors of leptons from Z decay and semi-leptonic top decay
c -- if flavors the same, then ambiguity over which comes from W and which from Z
      diffflav=.true.
      ST=0._dp

      if (runstring(4:5) == 'bt') then
         cut_id='btag'
         stop 'btag cuts not implemented for this process'
      elseif (runstring(4:5) == 'st') then
         cut_id='stee'
      endif
c -- write out cuts
      if (first) then
      write(*,*) '*----------- CMS FCNC search cuts -----------'
      write(*,*) '* '
      write(*,*) '*Lepton cuts:'
      write(*,*) '*pt,lepton > ', leptpt
      write(*,*) '* | eta,lepton| < ', leptrap
      write(*,*) mllmin, ' < m_ll < ',mllmax
      write(*,*) '*Lepton isolation: for all objects within ',Rleptisol,
     &     ' of a *lepton'
      write(*,*) '*sum (pt+E)/pt,lept <', Zisol, 'for lepton from Z'
      write(*,*) '*sum (pt+E)/pt,lept >', Wisol, 'for lepton from W'
      write(*,*) '*'
      write(*,*) '*Jet cuts'
      write(*,*) '*pt,jet >',ptjetmin
      write(*,*) '*| eta,j| <', etajetmax
      write(*,*) '* Jet isolation: R_jet,lept >', Rjlmin
      write(*,*) '*'
      write(*,*) '*Missing momentum pt,miss >', misspt
      write(*,*) '*'
      write(*,*) 'leptons have distinct flavors?',  diffflav
      if (cut_id == 'btag') then
      write(*,*) '*Top Mass cuts:'
      write(*,*) '*|m_Zj -mt| <', mZjcut
      write(*,*) '*|m_Wb-mt| <', mWbcut
      elseif (cut_id == 'stee') then
      write(*,*) mZjmin, ' < mZj < ', mZjmax
      write(*,*) mWbmin, ' < mWb < ', mWbmax
      write(*,*) 'ST > ', STmin

      endif
      write(*,*) '--------------------------------------------'
      first=.false.
      endif

c -- lepton index for t and tbar production
      if (nwz == +1) then
        leptindex=(/3,4,6/)
        nuindex=5
      elseif (nwz == -1) then
         leptindex=(/3,4,5/)
         nuindex=6
      endif

c     -- lepton cuts
      do j=1,3
         k=leptindex(j)
         ST=ST+pt(k,pjet)
c     -- pt and rap
         if ( (pt(k,pjet) <= leptpt) .or.
     &        (abs(etarap(k,pjet)) >= leptrap) ) then
            gencuts_WZjj=.true.
            return
         endif
      enddo
c -- dilepton mass cuts, for both l+l- pairings

c -- first establish which leptons reconstruct Z
      Zpair1=(/3,4/)
      Wpair1=(/5,6/)
      isol1=(/Zisol,Zisol,Wisol/)
      if (nwz == +1) then
         Zpair2=(/3,6/)
         Wpair2=(/4,5/)
         isol2=(/Zisol,Wisol,Zisol/)
      elseif (nwz == -1) then
         Zpair2=(/4,5/)
         Wpair2=(/3,6/)
         isol2=(/Wisol,Zisol,Zisol/)
      endif

      mll1=twomass(Zpair1(1),Zpair1(2),pjet)
      mll2=twomass(Zpair2(1),Zpair2(2),pjet)

      if (diffflav) then
c -- for different flavours, the W- and Z-pair are unambiguous
         Zpair=Zpair1
         Wpair=Wpair1
         isol=isol1
      if ( mll1 <= mllmin .or. mll1 >= mllmax ) then
         gencuts_WZjj=.true.
         return
      endif

      else
c -- for same flavours, check which l-l+ pair reconstructs Z mass best
         if ( abs(mll1-zmass) <= abs(mll2-zmass) ) then
            Zpair=Zpair1
            Wpair=Wpair1
            isol=isol1
            if ( mll1 <= mllmin .or. mll1 >= mllmax) then
               gencuts_WZjj=.true.
               return
            endif
         else
            Zpair=Zpair2
            Wpair=Wpair2
            isol=isol2
            if ( mll2 <= mllmin .or. mll2 >= mllmax) then
               gencuts_WZjj=.true.
               return
            endif
         endif
      endif


c---  identify the jets
      countjet=0
      do j=3,mxpart
         if ((plabel(j) == 'pp') .or. (plabel(j) == 'qj')
     &        .or. (plabel(j) =='bq') .or. (plabel(j) == 'ba')) then
            countjet=countjet+1
            jetindex(countjet)=j
         endif
      enddo

      incone=0._dp
      do j=1,3
c     --lepton isolation
         k=leptindex(j)
         do jj=1,3
            if (jj == j) cycle
            kk=leptindex(jj)
            if (R(pjet,k,kk) <= Rleptisol) then
               incone(j)=incone(j)+2._dp*pt(kk,pjet)
            endif
         enddo
         do jj=1,njets
            kk=jetindex(jj)
            if (R(pjet,k,kk) <= Rleptisol) then
               incone(j)=incone(j)+2._dp*pt(kk,pjet)
            endif
         enddo

      enddo

      if ( incone(1)/pt(leptindex(1),pjet) >= isol(1) .or.
     &     incone(2)/pt(leptindex(2),pjet) >= isol(2) .or.
     &     incone(3)/pt(leptindex(3),pjet) >= isol(3) ) then
         gencuts_WZjj=.true.
         return
      endif

c     --pt,miss
      ST=ST+pt(nuindex,pjet)
      if (pt(nuindex,pjet) <= misspt) then
         gencuts_WZjj=.true.
         return
      endif

c     -- jet cuts

c---  countjet will pick up the extra 'pp' needed for the real piece,
c---  therefore we should subtract 1 from this number
      if (countjet > njets) countjet=countjet-1

      if ((njets .ne. countjet) .and. (notag == 0)) then
         write(6,*) 'Something is wrong in gencuts.f -'
         write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
         stop
      endif

      if (njets <= 1) then
         gencuts_WZjj=.true.

         return
      endif

c -- jet isolation cuts
      do j=1,njets
         jj=jetindex(j)
         ST=ST+pt(jj,pjet)
         if ( (pt(jj,pjet) <= ptjetmin) .or.
     &        (abs(etarap(jj,pjet)) >= etajetmax) ) then
            gencuts_WZjj=.true.
            return
         endif

         do k=1,3
            kk=leptindex(k)
            if (R(pjet,jj,kk) <= Rjlmin) then
               gencuts_WZjj=.true.
               return
            endif
         enddo
      enddo

c -- reconstruct top mass using either jet and W/Z boson
      mZj1=threemass(Zpair(1),Zpair(2),jetindex(1),pjet)
      mZj2=threemass(Zpair(1),Zpair(2),jetindex(2),pjet)
      mWj1=threemass(Wpair(1),Wpair(2),jetindex(1),pjet)
      mWj2=threemass(Wpair(1),Wpair(2),jetindex(2),pjet)



      if (cut_id == 'stee') then
c -- "ST selection"

         if ( (mZj1 >= mZjmin .and. mZj1 <= mZjmax) .and.
     &        (mWj2 >= mWbmin .and. mWj2 <= mWbmax) ) then
         elseif ( (mZj2 >= mZjmin .and. mZj2 <= mZjmax) .and.
     &        (mWj1 >= mWbmin .and. mWj1 <= mWbmax) ) then
            continue
         else
            gencuts_WZjj=.true.
            return
         endif

       if ( ST <= STmin) then
          gencuts_WZjj=.true.
          return
       endif
      endif


       return
       end









