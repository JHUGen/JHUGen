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

      subroutine init_hist_Z
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

      do j=1,3
         if(j.eq.1) then
            prefix='p3'
         elseif(j.eq.2) then
            prefix='p4'
         elseif(j.eq.3) then
            prefix='p34'
         endif
         l=lastnbx(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.2d0,-6d0,6d0)
         call bookupeqbins(prefix(1:l)//'_eta',0.2d0,-6d0,6d0)
         call bookupeqbins(prefix(1:l)//'_pt' ,2d0,0d0,80d0)
         call bookupeqbins(prefix(1:l)//'_m'  ,5d0,0d0,100d0)
      enddo

      end




      subroutine analysis_Z(dsig0)
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
      real * 8 p3(4),p4(4),p34(4),y,eta,pt,mass
      integer j
      real * 8 prodvec2

      dsig  = dsig0

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         p3=phep(1:4,3)
         p4=phep(1:4,4)
         p34=p3+p4
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
         write(*,*) ' Now only NLO cross sections are implemented'
      endif

      call yetaptmassplot(p3,dsig,'p3')
      call yetaptmassplot(p4,dsig,'p4')
      call yetaptmassplot(p34,dsig,'p34')

      end


