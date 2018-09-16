      subroutine plots_stop_cfmt(p,wt,wt2)
      implicit none
      include 'types.f'
c--- routine to provide a plotting interface that produces the plots
c--- for both single top processes (2->2 and 2->3)
c--- of Campbell, Frederix, Maltoni and Tramontano
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      include 'nplot.f'
      include 'kprocess.f'
      include 'first.f'
      integer:: i,n,nplotmax,jet(2),jetswap,ibbar,inotb,ilight1,ilight2
      real(dp):: p(mxpart,4),wt,wt2,getet,
     & pttop,etatop,ytop,bwgt,
     & pt,etarap,yrap,ptthree,etarapthree,yrapthree,
     & wtbbar,wtnotb,wtlight1,wtlight2,
     & wtbbar2,wtnotb2,wtlight12,wtlight22,
     & ptmin,ptmax,ptbin,rapmin,rapmax,rapbin
      integer tag
      logical:: jetmerge
      common/jetmerge/jetmerge
      common/btagging/bwgt
      common/nplotmax/nplotmax
      integer, parameter:: tagbook=1, tagplot=2
!$omp threadprivate(/jetmerge/)
ccccc!$omp threadprivate(/nplotmax/)

c--- on the first call, initialize histograms
      if (first) then
        tag=tagbook
      pttop=-1._dp
      etatop=1d3
      ytop=1d3
      ibbar=6
      inotb=7
      ilight1=6
      ilight2=7
      goto 99
      else
        tag=tagplot
      endif
      
c--- check that this routine is being used for one of the envisaged
c--- processes and, if so, initialize variables accordingly      
      if     (kcase==kbq_tpq) then
        if (jets > 0) then
        do i=1,jets
          jet(i)=5+i ! jets start with parton 6 when removebr=.true.
        enddo
        pttop=ptthree(3,4,5,p)
        etatop=etarapthree(3,4,5,p)
        ytop=yrapthree(3,4,5,p)
      endif
      
      elseif (kcase==kqg_tbq) then
        if (jets > 0) then
        do i=1,jets
          jet(i)=4+i ! jets start with parton 5
        enddo
        pttop=pt(3,p)
        etatop=etarap(3,p)
        ytop=yrap(3,p)
      endif
      else
        write(6,*) 'This plotting routine is not suitable for'
      write(6,*) 'the process that you are calculating.'
      stop
      endif
      
c--- there are at most two jets; reorder them according to pt if two    
      if (jets == 2) then
        if (pt(jet(2),p) > pt(jet(1),p)) then
        jetswap=jet(1)
        jet(1)=jet(2)
        jet(2)=jetswap
      endif
      endif

c--- initialize all variabels to zero
      ibbar=0
      inotb=0
      ilight1=0
      ilight2=0

c--- Extract index to b~ (b for t~) in the 2->2 process
c--- Events should be plotted with weight = wt*bwgt instead of just wt
c--- (adapted from original code in nplotter.f by Z. Sullivan)
      if (kcase==kbq_tpq) then
c--- 1) Ascertain the identity of the jets
        if (jets > 0) then
c---     There are 2 cases: 
c---     If no merging, then only if jetlabel(i)=='pp' is it the b~.
          do i=1,jets
            if     (jetlabel(jet(i)-5)=='pp') then
            ibbar=jet(i)
            elseif (jetlabel(jet(i)-5)=='qj') then
            inotb=jet(i)
            endif        
          enddo
c---      If merging occurred, then it matters whether the merged
c---      jet was already a bq or not. 
          if (jetmerge) then
            ibbar=jet(1)   ! b~ was merged into jet(1)
          endif
      endif
c--- 2) Assign weights for histogramming
c---     there are three cases
        if     (jets == 1) then
          if     (jetlabel(1) == 'qj') then
c---     1) only one jet, that is definitely not a b
c---             jetlabel = (qj)
          ibbar=0
c           inotb already set above
          wtnotb=wt
          wtnotb2=wt2
        elseif (jetlabel(1) == 'pp') then
c---     2) only one jet, that could be a b
c---             jetlabel = (pp)
c           ibbar already set above
          inotb=ibbar
          wtbbar=wt*bwgt
          wtnotb=wt*(1._dp-bwgt)
          wtbbar2=wt2*bwgt**2
          wtnotb2=wt2*(1._dp-bwgt)**2
          endif        
c---     3) two jets, one of which could be a b
c---             jetlabel = (qj,pp) OR (pp,qj)
        elseif (jets == 2) then
c         ibbar already set above
c         inotb already set above
          wtbbar=wt*bwgt
        wtnotb=wt*bwgt
        wtbbar2=wt2*bwgt**2
        wtnotb2=wt2*bwgt**2
          ilight1=jet(1) ! they are ordered according to pt
          ilight2=jet(2) ! they are ordered according to pt
        wtlight1=wt*(1._dp-bwgt)
        wtlight2=wtlight1
        wtlight12=wt2*(1._dp-bwgt)**2
        wtlight22=wtlight12
        elseif (jets == 0) then
c---   nothing to set: all variables=0, only occurs when notag=1
          continue
        else
          write(6,*) 'Error: there should be 1 or 2 jets, instead ',jets
        stop
        endif
      endif

c--- Simpler manipulations for the 2->3 process
      if (kcase==kqg_tbq) then
        ibbar=4
      wtbbar=wt
      wtbbar2=wt2
      if (jets > 0) then
        ilight1=jet(1)
        wtlight1=wt
        wtlight12=wt2
      endif
      if (jets > 1) then
        ilight2=jet(2)
        wtlight2=wt
        wtlight22=wt2
      endif
c--- inotb should remain equal to zero              
      endif

   99 continue
      n=1

c--- fill plots
c--- available veriables are
c---   pttop,eta,ytop: pt, pseudo-rapidity and rapidity of the top quark
      
      call bookplot(n,tag,'pt(top) ',pttop,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1   
      call bookplot(n,tag,'eta(top)',etatop,wt,wt2,-8._dp,8._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y(top)  ',ytop,wt,wt2,-8._dp,8._dp,0.5_dp,'lin')
      n=n+1
      
c--- more complicated procedure for these plots, due to 2->2 process
c---  pt and rapidity of the b-quark and other jets 
c---  Note that the order is important (due to histogram finalizing):
c---     all plots for jet 1, all for bottom quark, all for jet 2     

c--- set the ranges and bin sizes here
      ptmin=0._dp
      ptmax=200._dp
      ptbin=2._dp
      rapmin=-5._dp
      rapmax=+5._dp
      rapbin=0.2_dp      

c--- 1) leading jet that isn't a b
c---    first account for events with only one non b jet (always occurs)
c---    (always occurs, except if notag=1)
c--- [pt]
      if (inotb .ne. 0) then
      call bookplot(n,tag,'pt(jet 1)',
     & pt(inotb,p),wtnotb,wtnotb2,ptmin,ptmax,ptbin,'log')
      endif
c---    then account for events with two non b jets (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'pt(jet 1)',
     & pt(ilight1,p),wtlight1,wtlight12,ptmin,ptmax,ptbin,'log')
      endif
      n=n+1
c--- [rapidity]
      if (inotb .ne. 0) then
      call bookplot(n,tag,'eta(jet 1)',
     & etarap(inotb,p),wtnotb,wtnotb2,rapmin,rapmax,rapbin,'lin')
      endif
c---    then account for events with two non b jets (only if ilight1>0)
      if (ilight1 .ne. 0) then
      call bookplot(n,tag,'eta(jet 1)',
     & etarap(ilight1,p),wtlight1,wtlight12,rapmin,rapmax,rapbin,'lin')
      endif
      n=n+1

c--- 2) b jet (only if ibbar>0)
c--- [pt]
      if (ibbar .ne. 0) then
      call bookplot(n,tag,'pt(bottom)',
     & pt(ibbar,p),wtbbar,wtbbar2,ptmin,ptmax,ptbin,'log')
      endif
      n=n+1
c--- [pt] with |rapidity| < 2.8  (plot for CDF)
      if (ibbar .ne. 0) then
        if ((abs(etarap(ibbar,p)) < 2.8_dp).or.(tag == tagbook)) then
        call bookplot(n,tag,'pt(bottom) with |eta(bottom)| < 2.8',
     &   pt(ibbar,p),wtbbar,wtbbar2,ptmin,ptmax,ptbin,'log')
        endif
      endif
      n=n+1
c--- [rapidity]
      if (ibbar .ne. 0) then
      call bookplot(n,tag,'eta(bottom)',
     & etarap(ibbar,p),wtbbar,wtbbar2,rapmin,rapmax,rapbin,'lin')
      endif
      n=n+1
c--- [rapidity] < 2.8 and pt > 20 GeV (plot for acceptance in CDF)
      if (ibbar .ne. 0) then
        if ((pt(ibbar,p) > 20._dp).or.(tag == tagbook)) then
        call bookplot(n,tag,'eta(bottom) with pt(bottom)>20 GeV',
     &   etarap(ibbar,p),wtbbar,wtbbar2,-2.8_dp,2.801_dp,rapbin,'lin')
        endif
      endif
      n=n+1
c--- [rapidity] < 2.8 and pt > 20 GeV (plot for acceptance in CDF)
      if (ibbar .ne. 0) then
        if ((getet(p(ibbar,4),p(ibbar,1),p(ibbar,2),p(ibbar,3))>20._dp)
     &      .or.(tag == tagbook)) then
        call bookplot(n,tag,'eta(bottom) with Et(bottom)>20 GeV',
     &   etarap(ibbar,p),wtbbar,wtbbar2,-2.8_dp,2.801_dp,rapbin,'lin')
        endif
      endif
      n=n+1

c--- 3) subleading jet that isn't a b (only if ilight2>0)
c--- [pt]
      if (ilight2 .ne. 0) then
      call bookplot(n,tag,'pt(jet 2)',
     & pt(ilight2,p),wtlight2,wtlight22,ptmin,ptmax,ptbin,'log')
      endif
      n=n+1
c--- [rapidity]
      if (ilight2 .ne. 0) then
      call bookplot(n,tag,'eta(jet 2)',
     & etarap(ilight2,p),wtlight2,wtlight22,rapmin,rapmax,rapbin,'lin')
      endif
      n=n+1

      
c--- copied from nplotter.f      
      n=n-1

      if (n > maxhisto) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > ',maxhisto,', which is the built-in maximum.'
      write(6,*) 'To use more histograms, change the value of the'
      write(6,*) 'constant MAXHISTO in src/Inc/nplot.f then do:'
      write(6,*)
      write(6,*) ' make clean; make        to recompile from scratch.'
        write(6,*)
      stop
      endif

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
            
      return
      end
      
