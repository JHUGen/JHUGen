C ----------------------------------------------------------------- C
C - This is a parton level only analysis for ttbar production     - C
C - tops bs Ws etc are constructed using MC truth only - no jets  - C
C - The MC truth reconstruction has been tested (see sanity check - C
C - code in the analysis below).                                  - C
C - Since it is parton level you need to comment out from         - C
C - CALL HWDHOB down to CALL HWDHOB inclusive in main-HERWIG.f .  - C
C - Also at some point, the showering went into what looked like  - C
C - an infinite loop after 137K events - gdb said it was in       - C
C - HWHGUP. The same glitch did not occur with *** herwig6520.f *** C
C - I also eliminated the analysis as a possible cause (it occurs - C
C - with HWANAL removed). I did not see anything fishy with the   - C
C - Tevatron, semileptonic event that got caught.                 - C
C ----------------------------------------------------------------- C

c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist_KN
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer lastnbx,j,l
      character * 20 prefix
      integer nbins
      parameter (nbins=11)
      real * 8 pT_tt_bins(nbins+1)
      data pT_tt_bins/  0d0, 10d0, 25d0, 50d0,100d0,
     1     150d0,200d0,250d0,300d0,400d0,600d0, 900d0/          
      real * 8 m_tt_bins(nbins+1)
      data m_tt_bins/ 320d0,360d0,400d0,450d0,500d0,
     1     550d0,600d0,650d0,700d0,800d0,900d0,1000d0/          
      external lastnbx

      call inihists

      do j=1,18
         if(j.eq.1) then
            prefix='t'
         elseif(j.eq.2) then
            prefix='tb'
         elseif(j.eq.3) then
            prefix='t-pt5'
         elseif(j.eq.4) then
            prefix='tb-pt5'
         elseif(j.eq.5) then
            prefix='t-pt20'
         elseif(j.eq.6) then
            prefix='tb-pt20'
         elseif(j.eq.7) then
            prefix='btop'
         elseif(j.eq.8) then
            prefix='btop-pt5'
         elseif(j.eq.9) then
            prefix='btop-pt20'
         elseif(j.eq.10) then
            prefix='bbtop'
         elseif(j.eq.11) then
            prefix='bbtop-pt5'
         elseif(j.eq.12) then
            prefix='bbtop-pt20'
         elseif(j.eq.13) then
            prefix='lwp'
         elseif(j.eq.14) then
            prefix='lwp-pt5'
         elseif(j.eq.15) then
            prefix='lwp-pt20'
         elseif(j.eq.16) then
            prefix='lwm'
         elseif(j.eq.17) then
            prefix='lwm-pt5'
         elseif(j.eq.18) then
            prefix='lwm-pt20'
         endif
         l=lastnbx(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_eta',0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_pt' ,5d0,0d0,400d0)
         call bookupeqbins(prefix(1:l)//'_m'  ,2d0,-0.5d0,201.5d0)
      enddo

      end




      subroutine analysis_KN(dsig0)
      implicit none
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      integer jpref
      character * 20 prefix(18)
      common/ccccprefix/jpref,prefix
      real * 8  dsig0,dsig
      logical   ini
      data      ini/.true./
      save      ini
      integer   ihep                ! HEPEVT index.
      real * 8 p_top(4),p_tb(4),p_wp(4),p_wm(4),p_lwp(4),p_lwm(4),
     1         p_b(4),p_bb(4),y,eta,pt,mass,ptt
      integer j
      real * 8 prodvec2

      dsig  = dsig0

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         p_top=phep(1:4,3)
         p_tb=phep(1:4,4)
         p_wp=phep(1:4,5)
         p_wm=phep(1:4,6)
         p_lwp=phep(1:4,8)
         p_lwm=phep(1:4,9)
         p_b=phep(1:4,11)
         p_bb=phep(1:4,12)
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
         write(*,*) ' Now only NLO cross sections are implemented'
      endif

      if(nhep.eq.13) then
         ptt=phep(1,13)**2+phep(2,13)**2
         do j=7,12
            ptt=min(ptt,prodvec2(phep(:,j),phep(:,13)))
         enddo
         ptt=sqrt(ptt)
      else
         ptt=0
      endif

      call yetaptmassplot(p_top,dsig,'t')
      if(ptt.gt.5) call yetaptmassplot(p_top,dsig,'t-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_top,dsig,'t-pt20')

      call yetaptmassplot(p_tb,dsig,'tb')
      if(ptt.gt.5) call yetaptmassplot(p_tb,dsig,'tb-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_tb,dsig,'tb-pt20')

      call yetaptmassplot(p_b,dsig,'btop')
      if(ptt.gt.5) call yetaptmassplot(p_b,dsig,'btop-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_b,dsig,'btop-pt20')

      call yetaptmassplot(p_bb,dsig,'bbtop')
      if(ptt.gt.5) call yetaptmassplot(p_bb,dsig,'bbtop-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_bb,dsig,'bbtop-pt20')

      call yetaptmassplot(p_lwp,dsig,'lwp')
      if(ptt.gt.5) call yetaptmassplot(p_lwp,dsig,'lwp-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_lwp,dsig,'lwp-pt20')

      call yetaptmassplot(p_lwm,dsig,'lwm')
      if(ptt.gt.5) call yetaptmassplot(p_lwm,dsig,'lwm-pt5')
      if(ptt.gt.20) call yetaptmassplot(p_lwm,dsig,'lwm-pt20')

      end


C ----------------------------------------------- C
C - SANITY CHECK OF TRUTH SHOWER RECONSTRUCTION - C
C ----------------------------------------------- C
C - p_tmp is the total momentum in minus the total momentum out.
C - If the sum of the absolute values of p_tmp are greater than
C - epsilon (1 GeV) a warning plus debugging output is printed
C - on the screen (see below). At Tevatron with 300K events I
C - see in Herwig messages displayed for epsilon equal to 5 GeV
C - every 30K events - n.b. these violations can be directly
C - and totally attributed exclusively to differences between
C - the Wminus momentum reconstructed and it's true value (see
C - below), what's more 90% of them are associated to the W
C - decaying to strange+charm which then shower - the strange
C - and charm momenta are the things that change quite a bit
C - between W decay and post-showering. Pythia is totally fine
C - for 300K events - no messages. Possibly this is a HW issue then.

      subroutine shower_sanity_checks(p_top,p_tba,p_beam,p_Wp,p_Wm,
     1                                p_b  ,p_bba,p_jet)
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      real * 8  p_top(4),p_tba(4),p_beam(4),p_tmp(4)
      real * 8  p_Wp(4) ,p_Wm(4) ,p_b(4)   ,p_bba(4),p_jet(4)
c$$$      real * 8           y_top,eta_top,pt_top,m_top,
c$$$     1                   y_tba,eta_tba,pt_tba,m_tba
c$$$      real * 8           y_jet,eta_jet,pt_jet,m_jet ! Powheg _radiated_ jet 
      real * 8  tmp
      integer   n_sanity_violations,n_event
      data      n_sanity_violations,n_event/0,0/
      save      n_sanity_violations,n_event
      integer   ixx

c$$$      call getyetaptmass(p_top,y_top,eta_top,pt_top,m_top)
c$$$      call getyetaptmass(p_tba,y_tba,eta_tba,pt_tba,m_tba)
c$$$      call getyetaptmass(p_jet,y_jet,eta_jet,pt_jet,m_jet)

      n_event=n_event+1

      p_tmp=p_top+p_tba+p_jet+p_beam-phep(1:4,1)-phep(1:4,2)
      tmp=0d0
      do ixx=1,4
         tmp=tmp+abs(p_tmp(ixx))
      enddo
      if(abs(tmp).gt.5d0) then
         n_sanity_violations=n_sanity_violations+1
         if(n_sanity_violations.le.10) then
            write(*,*) ''
            write(*,*) ''
            if(whcprg.eq.'HERWIG') then
               write(*,*) 'HERWIG analysis'
CC                  call HWUEPR
CC                  call HWUBPR
            else if(whcprg.eq.'PYTHIA') then
               write(*,*) 'PYTHIA analysis'
CC                  call pylist(5)
            else if(whcprg.eq.'LHE') then
               write(*,*) 'LHE analysis'
               do ixx=1,nhep
                  write(*,*) 'ihep =',ixx,' id = ',idhep(ixx),
     1                       ' mom. = ',phep(1:4,ixx)
               enddo
            endif
            write(*,*) ''
            write(*,*) 'Consistency err: momentum in != momentum out'
            write(*,*) 'Event number ',n_event
            write(*,*) 'Debugging output follows.'
            write(*,*) 'nup   = ',nup
            write(*,*) 'p_in  = ',phep(1:4,1)+phep(1:4,2)
            write(*,*) 'p_out = ',p_top+p_tba+p_jet+p_beam
            
            if(whcprg.eq.'PYTHIA') then
               p_tmp=phep(1:4,7)
               write(*,*) 'p_top true = ',p_tmp
            else if(whcprg.eq.'HERWIG') then
               do ixx=1,nhep
                  if(isthep(ixx).eq.155.and.idhep(ixx).eq.6) then
                     p_tmp=phep(1:4,ixx)
                  endif
               enddo
               write(*,*) 'p_top true = ',p_tmp
            endif
            write(*,*) 'p_top      = ',p_top
c$$$            write(*,*) 'm_top      = ',m_top
            
            do ixx=1,nhep
               if(whcprg.eq.'PYTHIA') then
                  if(isthep(ixx).eq.2.and.idhep(ixx).eq.24) then
                     p_tmp=phep(1:4,ixx)
                     write(*,*) 'p_Wp true  = ',p_tmp
                  endif
               else if(whcprg.eq.'HERWIG') then
                  if(isthep(ixx).eq.155.and.idhep(ixx).eq.24) then
                     p_tmp=phep(1:4,ixx)
                     write(*,*) 'p_Wp true  = ',p_tmp
                  endif
               endif
            enddo
            write(*,*) 'p_Wp       = ',p_Wp
            
            if(whcprg.eq.'PYTHIA') then
               p_tmp=phep(1:4,11)
               write(*,*) 'p_b true   = ',p_tmp
            else if(whcprg.eq.'HERWIG') then
               do ixx=1,nhep
                  if(isthep(ixx).eq.124.and.idhep(ixx).eq.5) then
                     p_tmp=phep(1:4,ixx)
                  endif
               enddo
               write(*,*) 'p_b true   = ',p_tmp
            endif
            write(*,*) 'p_b        = ',p_b
            
            
            if(whcprg.eq.'PYTHIA') then
               p_tmp=phep(1:4,8)
               write(*,*) 'p_tba true = ',p_tmp
            else if(whcprg.eq.'HERWIG') then
               do ixx=1,nhep
                  if(isthep(ixx).eq.155.and.idhep(ixx).eq.-6) then
                     p_tmp=phep(1:4,ixx)
                  endif
               enddo
               write(*,*) 'p_tba true = ',p_tmp
            endif
            write(*,*) 'p_tba      = ',p_tba
c$$$            write(*,*) 'm_tba      = ',m_tba
            
            do ixx=1,nhep
               if(whcprg.eq.'PYTHIA') then
                  if(isthep(ixx).eq.2.and.idhep(ixx).eq.-24) then
                     p_tmp=phep(1:4,ixx)
                     write(*,*) 'p_Wm true  = ',p_tmp
                  endif
               else if(whcprg.eq.'HERWIG') then
                  if(isthep(ixx).eq.155.and.idhep(ixx).eq.-24) then
                     p_tmp=phep(1:4,ixx)
                     write(*,*) 'p_Wm true  = ',p_tmp
                  endif
               endif
            enddo
            write(*,*) 'p_Wm       = ',p_Wm
            
            if(whcprg.eq.'PYTHIA') then
               p_tmp=phep(1:4,13)
               write(*,*) 'p_bba true = ',p_tmp
            else if(whcprg.eq.'HERWIG') then
               do ixx=1,nhep
                  if(isthep(ixx).eq.124.and.idhep(ixx).eq.-5) then
                     p_tmp=phep(1:4,ixx)
                  endif
               enddo
               write(*,*) 'p_bba true = ',p_tmp
            endif
            write(*,*) 'p_bba      = ',p_bba
            
            if(whcprg.eq.'PYTHIA') then
               write(*,*) 'p_jet true = ',phep(1:4,9)
            else if(whcprg.eq.'HERWIG') then
               write(*,*) 'p_jet true = ',phep(1:4,jdahep(1,9))
            endif
            write(*,*) 'p_jet      = ',p_jet
            
            write(*,*) 'p_beam      = ',p_beam

         else if(n_sanity_violations.eq.11) then
            write(*,*) 'Analysis found 10 momentum conservation'
            write(*,*) 'inconsistencies - no more debug ouput'
            write(*,*) 'will be given on these, just warnings.'
         else if(mod(n_sanity_violations,20).eq.0) then
            write(*,*) 'Analysis has now found',
     1           n_sanity_violations,' momentum inconsistencies.'
         endif
      endif

      end



