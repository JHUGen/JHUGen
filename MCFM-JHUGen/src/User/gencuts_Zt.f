      logical function gencuts_Zt(pjet,njets)
c -- these cuts taken from CMS FCNC search, arxiv:1208.0957
c -- R. Rontsch 2103-03-01
      implicit none
      include 'leptcuts.f'
      include 'jetcuts.f'
      include 'constants.f'
      include 'plabel.f'
      include 'jetlabel.f'
      include 'masses.f'
      include 'nwz.f'
      include 'runstring.f'
      include 'notag.f'
      integer njets
      double precision pjet(mxpart,4)
      integer leptindex(3), jetindex(mxpart),countjet
      integer j,k,jj,kk, Zpair1(2),Zpair2(2),Wpair1(2),Wpair2(2),
     &     Wpair(2),Zpair(2),nuindex
      logical first, diffflav
      double precision mll1,mll2,etarap,twomass,threemass,pt,r
      double precision mllmax,mllmin,Rleptisol,Zisol,Wisol,
     &     mZjcut,mWbcut,mZj,mZjold,mWb,incone(3),mWbmin,mWbmax,
     &     mZjmin,mZjmax,STmin,ST,isol(3),isol1(3),isol2(3)
      parameter (mllmin=60d0,mllmax=120d0,Rleptisol=0.3d0,Zisol=0.125d0,
     &     Wisol=0.1d0,mZjcut=25d0,mWbcut=35d0,mZjmin=100d0,
     &     mZjmax=250d0,mWbmin=100d0,mWbmax=250d0,STmin=250d0)
      character*4 cut_id
      data first/.true./

      gencuts_Zt=.false.
c -- allows for same flavors of leptons from Z decay and semi-leptonic top decay
c -- if flavors the same, then ambiguity over which comes from W and which from Z
      diffflav=.true.
      ST=0d0

c -- CMS "S_T" cuts and "b-tag" cuts
      if (runstring(4:5) .eq. 'bt') then
         cut_id='btag'
      elseif (runstring(4:5) .eq. 'st') then
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
      if (cut_id .eq. 'btag') then
      write(*,*) '*Top Mass cuts:'
      write(*,*) '*|m_Zj -mt| <', mZjcut
      write(*,*) '*|m_Wb-mt| <', mWbcut
      elseif (cut_id .eq. 'stee') then
      write(*,*) mZjmin, ' < mZj < ', mZjmax
      write(*,*) mWbmin, ' < mWb < ', mWbmax
      write(*,*) 'ST > ', STmin

      endif
      write(*,*) '--------------------------------------------'
      first=.false.
      endif

c -- lepton index for t and tbar production
      if (nwz .eq. +1) then
        leptindex=(/3,4,6/)
        nuindex=5
      elseif (nwz .eq. -1) then
         leptindex=(/3,4,5/)
         nuindex=6
      endif

c     -- lepton cuts
      do j=1,3
         k=leptindex(j)
         ST=ST+pt(k,pjet)
c     -- pt and rap
         if ( (pt(k,pjet) .le. leptpt) .or.
     &        (abs(etarap(k,pjet)) .ge. leptrap) ) then
            gencuts_Zt=.true.
            return
         endif
      enddo

c -- dilepton mass cuts, for both l+l- pairings
c -- first establish which leptons reconstruct Z
      Zpair1=(/3,4/)
      Wpair1=(/5,6/)
      isol1=(/Zisol,Zisol,Wisol/)
      if (nwz .eq. +1) then
         Zpair2=(/3,6/)
         Wpair2=(/4,5/)
         isol2=(/Zisol,Wisol,Zisol/)
      elseif (nwz .eq. -1) then
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
      if ( mll1 .le. mllmin .or. mll1 .ge. mllmax ) then
         gencuts_Zt=.true.
         return
      endif

      else
c -- for same flavours, check which l-l+ pair reconstructs Z mass best
         if ( abs(mll1-zmass) .le. abs(mll2-zmass) ) then
            Zpair=Zpair1
            Wpair=Wpair1
            isol=isol1
            if ( mll1 .le. mllmin .or. mll1 .ge. mllmax) then
               gencuts_Zt=.true.
               return
            endif
         else
            Zpair=Zpair2
            Wpair=Wpair2
            isol=isol2
            if ( mll2 .le. mllmin .or. mll2 .ge. mllmax) then
               gencuts_Zt=.true.
               return
            endif
         endif
      endif


c---  identify the jets
      countjet=0
      do j=3,mxpart
         if ((plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     &        .or. (plabel(j) .eq.'bq') .or. (plabel(j) .eq. 'ba')) then
            countjet=countjet+1
            jetindex(countjet)=j
         endif
      enddo

      incone=0d0
      do j=1,3
c     -- lepton isolation
         k=leptindex(j)
         do jj=1,3
            if (jj .eq. j) cycle
            kk=leptindex(jj)
            if (R(pjet,k,kk) .le. Rleptisol) then
               incone(j)=incone(j)+2d0*pt(kk,pjet)
            endif
         enddo
         do jj=1,njets
            kk=jetindex(jj)
            if (R(pjet,k,kk) .le. Rleptisol) then
               incone(j)=incone(j)+2d0*pt(kk,pjet)
            endif
         enddo

      enddo

      if ( incone(1)/pt(leptindex(1),pjet) .ge. isol(1) .or.
     &     incone(2)/pt(leptindex(2),pjet) .ge. isol(2) .or.
     &     incone(3)/pt(leptindex(3),pjet) .ge. isol(3) ) then
         gencuts_Zt=.true.
         return
      endif

c     --pt,miss
      ST=ST+pt(nuindex,pjet)
      if (pt(nuindex,pjet) .le. misspt) then
         gencuts_Zt=.true.
         return
      endif


c     -- jet cuts

c---  countjet will pick up the extra 'pp' needed for the real piece,
c---  therefore we should subtract 1 from this number
      if (countjet .gt. njets) countjet=countjet-1

      if ((njets .ne. countjet) .and. (notag .eq. 0)) then
         write(6,*) 'Something is wrong in gencuts.f -'
         write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
         stop
      endif

      if (njets .le. 1) then
         gencuts_Zt=.true.
         return
      endif

c -- jet isolation cuts
      do j=1,njets
         jj=jetindex(j)
         ST=ST+pt(jj,pjet)
         if ( (pt(jj,pjet) .le. ptjetmin) .or.
     &        (abs(etarap(jj,pjet)) .ge. etajetmax) ) then
            gencuts_Zt=.true.
            return
         endif

         do k=1,3
            kk=leptindex(k)
            if (R(pjet,jj,kk) .le. Rjlmin) then
               gencuts_Zt=.true.
               return
            endif
         enddo
      enddo

c     -- find jet that best reconstructs top mass
      mZjold=1d5
      do j=1,njets
         jj=jetindex(j)
         if (jetlabel(j) .eq. 'bq' .or. jetlabel(j) .eq.'ba') cycle
         mZj=threemass(Zpair(1),Zpair(2),jj,pjet)
         if ( abs(mZj-mt) .le. abs(mZjold-mt)) then
            mZjold=mZj
         endif
      enddo

      do j=1,njets
         jj=jetindex(j)
         if (jetlabel(j) .eq. 'bq' .or. jetlabel(j) .eq.'ba') then
            mWb=threemass(Wpair(1),Wpair(2),jj,pjet)
         endif
      enddo


c     -- "b-tag selection"
      if (cut_id .eq. 'btag') then

c     -- cut on mZj
      if ( abs(mZjold-mt) .ge. mZjcut) then
         gencuts_Zt=.true.
         return
      endif

c     -- cut on mWb
      if ( abs(mWb-mt) .ge. mWbcut) then
         gencuts_Zt=.true.
         return
      endif


      elseif (cut_id .eq. 'stee') then


c -- "ST selection"

      if (mZjold .le. mZjmin .or. mZjold .ge. mZjmax) then
          gencuts_Zt=.true.
          return
       endif

       if (mWb .le. mWbmin .or. mWb .ge. mWbmax) then
          gencuts_Zt=.true.
          return
       endif

       if ( ST .le. STmin) then
          gencuts_Zt=.true.
          return
       endif
      endif

       return
       end









