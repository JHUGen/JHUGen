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

      subroutine init_hist_ST_sch_dk
      implicit none
      include 'types.f'
      
      include  'LesHouches.h'
      include 'pwhg_math.h'
      integer:: lastnbx,j,l
      character * 20 prefix
      integer:: nbins
      parameter (nbins=11)
      real * 8 pT_tt_bins(nbins+1)
      data pT_tt_bins/  0._dp, 10._dp, 25._dp, 50._dp,100._dp,
     1     150._dp,200._dp,250._dp,300._dp,400._dp,600._dp, 900._dp/          
      real * 8 m_tt_bins(nbins+1)
      data m_tt_bins/ 320._dp,360._dp,400._dp,450._dp,500._dp,
     1     550._dp,600._dp,650._dp,700._dp,800._dp,900._dp,1000._dp/          
      external lastnbx

      call inihists

      do j=1,6
         if(j==1) then
            prefix='t'
         elseif(j==2) then
            prefix='b'
         elseif(j==3) then
            prefix='W'
         elseif(j==4) then
            prefix='l'
         elseif(j==5) then
            prefix='a'
         elseif(j==6) then
            prefix='btop'
          endif
         l=lastnbx(prefix)
         call bookupeqbins(prefix(1:l)//'_y'  ,0.2_dp,-4._dp,4._dp)
         call bookupeqbins(prefix(1:l)//'_eta',0.2_dp,-4._dp,4._dp)
         call bookupeqbins(prefix(1:l)//'_pt' ,5._dp,0._dp,400._dp)
         call bookupeqbins(prefix(1:l)//'_m'  ,2._dp,-0.5_dp,201.5_dp)
      enddo

      end




      subroutine analysis_ST_sch_dk(dsig0)
      implicit none
      include 'types.f'
      
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      integer:: jpref
      character * 20 prefix(18)
      common/ccccprefix/jpref,prefix
      real * 8  dsig0,dsig
      logical::   ini
      data      ini/.true./
      save      ini
      integer::   ihep                ! HEPEVT index.
      real * 8 p_top(4),p_b(4),p_w(4),p_l(4),p_a(4),p_btop(4),
     1         y,eta,pt,mass,ptb
      integer:: j
      real * 8 prodvec2

      dsig  = dsig0

      if(whcprg=='NLO'.or.whcprg=='LHE') then
         p_top=phep(1:4,3)
         p_b=phep(1:4,8)
         p_w=phep(1:4,4)
         p_l=phep(1:4,5)
         p_a=phep(1:4,6)
         p_btop=phep(1:4,7)
C --------------------------------------------- C
C - LHE PARTICLE TOP RECONSTRUCTION: MC TRUTH - C
C --------------------------------------------- C
      else
         write(*,*) ' Now only NLO cross sections are implemented'
      endif

      if(nhep==9) then
         ptb=phep(1,9)**2+phep(2,9)**2
         do j=5,8
            ptb=min(ptb,prodvec2(phep(:,j),phep(:,9)))
         enddo
         ptb=sqrt(ptb)
      else
         ptb=0
      endif

      call yetaptmassplot(p_top,dsig,'t')

      call yetaptmassplot(p_b,dsig,'b')

      call yetaptmassplot(p_w,dsig,'W')

      call yetaptmassplot(p_l,dsig,'l')

      call yetaptmassplot(p_a,dsig,'a')

      call yetaptmassplot(p_btop,dsig,'btop')

      end


