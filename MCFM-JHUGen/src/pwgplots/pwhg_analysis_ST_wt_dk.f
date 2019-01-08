C ----------------------------------------------------------------- C
C - This is a parton level only analysis for s-channel single top - C
C - production including the top quark decay                      - C
C - tops bs Ws etc are constructed using MC truth only - no jets  - C
C - The MC truth reconstruction has not been tested               - C
C ----------------------------------------------------------------- C

c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist_ST_wt_dk
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

      do j=1,8
         if(j.eq.1) then
            prefix='t'
         elseif(j.eq.2) then
            prefix='W'
         elseif(j.eq.3) then
            prefix='Wtop'
         elseif(j.eq.4) then
            prefix='ltop'
         elseif(j.eq.5) then
            prefix='atop'
         elseif(j.eq.6) then
            prefix='btop'
         elseif(j.eq.7) then
            prefix='l'
         elseif(j.eq.8) then
            prefix='a'
          endif
         l=lastnbx(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_eta',0.2d0,-4d0,4d0)
         call bookupeqbins(prefix(1:l)//'_pt' ,5d0,0d0,400d0)
         call bookupeqbins(prefix(1:l)//'_m'  ,2d0,-0.5d0,201.5d0)
      enddo

      end




      subroutine analysis_ST_wt_dk(dsig0)
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
      real * 8 p_top(4),p_wtop(4),p_ltop(4),p_atop(4),p_btop(4),
     &   p_w(4),p_l(4),p_a(4),y,eta,pt,mass,ptWt
      integer j
      real * 8 prodvec2

      dsig  = dsig0

      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then
         p_top=phep(1:4,3)
         p_wtop=phep(1:4,4)
         p_w=phep(1:4,5)
         p_l=phep(1:4,6)
         p_a=phep(1:4,7)
         p_ltop=phep(1:4,8)
         p_atop=phep(1:4,9)
         p_btop=phep(1:4,10)
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
         write(*,*) ' Now only NLO cross sections are implemented'
      endif

      if(nhep.eq.11) then
         ptWt=phep(1,11)**2+phep(2,11)**2
c         do j=5,8
c            ptb=min(ptb,prodvec2(phep(:,j),phep(:,11)))
c         enddo
         ptWt=sqrt(ptWt)
      else
         ptWt=0
      endif

      call yetaptmassplot(p_top,dsig,'t')

      call yetaptmassplot(p_w,dsig,'W')

      call yetaptmassplot(p_wtop,dsig,'Wtop')

      call yetaptmassplot(p_ltop,dsig,'ltop')

      call yetaptmassplot(p_atop,dsig,'atop')

      call yetaptmassplot(p_btop,dsig,'btop')

      call yetaptmassplot(p_l,dsig,'l')

      call yetaptmassplot(p_a,dsig,'a')

      end


