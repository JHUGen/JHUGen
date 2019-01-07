      subroutine nplotter_ttbar(p,wt,wt2,switch)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'masses.f'
      include 'plabel.f'
      include 'process.f'
      include 'outputflags.f'
      double precision p(mxpart,4),wt,wt2,yrap,pt,ptjet(mxpart),
     & pttwo,mll,ptl,yl,ptttb,yttb,dot,tiny,
     & ptt,pttb,yt,ytb,mttb,mwp,mwm,mlb
      double precision plep(4),ptop(4),pleprest(4),ptoprest(4),mj1j2
      integer switch,n,nplotmax,j,
     & jetindex(mxpart),iorder(mxpart),ijet
      character*4 tag
      logical dilepton,failed
      parameter(tiny=1d-8)
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,mj1j2)
      
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- determine whether or not process represents dilepton channel
      if ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_ldk')
     &.or.(case .eq. 'tt_bbu') .or. (case .eq. 'tt_udk') ) then
        dilepton=.true.
      else
        dilepton=.false.
      endif
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        mj1j2=0d0
        goto 99
      else
c--- Add event in histograms
        tag='plot'
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

       ijet=0
       do j=3,9
         if ((plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'bq')
     &   .or.(plabel(j) .eq. 'ba')) then
           if (p(j,4) .gt. tiny) then
             ijet=ijet+1
             jetindex(ijet)=j
             ptjet(ijet)=pt(j,p)
           endif
         endif
       enddo   
       call arraysort(ijet,ptjet,iorder)

       if (ijet .ge. 2) then
c--- Two hardest jets are now indexed by
c--- jetindex(iorder(1)) and jetindex(iorder(2))
       mj1j2=((p(jetindex(iorder(1)),4)+p(jetindex(iorder(2)),4))**2
     &       -(p(jetindex(iorder(1)),1)+p(jetindex(iorder(2)),1))**2
     &       -(p(jetindex(iorder(1)),2)+p(jetindex(iorder(2)),2))**2
     &       -(p(jetindex(iorder(1)),3)+p(jetindex(iorder(2)),3))**2)
       mj1j2=dsqrt(max(mj1j2,0d0))
       else
       mj1j2=-1d0
       endif
       
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/dfloat(itmx))  
      endif

c--- "n" will count the number of histograms
      n=nextnplot              

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

c--- This section tailored for comparison with CMS PAS TOP-11-013

      if (dilepton) then
c--- dilepton=specific plots
                   
c--- pt(l+ and l-)
      call bookplot(n,tag,'pt(l+ and l-)',pt(4,p),wt,wt2,
     & 0d0,400d0,10d0,'log')
      call bookplot(n,tag,'pt(l+ and l-)',pt(7,p),wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1

c--- eta(l+ and l-)
      call bookplot(n,tag,'eta(l+ and l-)',yrap(4,p),wt,wt2,
     & -2.4d0,2.4d0,0.2d0,'lin')
      call bookplot(n,tag,'eta(l+ and l-)',yrap(7,p),wt,wt2,
     & -2.4d0,2.4d0,0.2d0,'lin')
      n=n+1

c--- pt(l+,l-)
      call bookplot(n,tag,'pt(l+,l-)',pttwo(4,7,p),wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1

c--- m(l+,l-)
      mll=dsqrt(max(2d0*dot(p,4,7),0d0))
      call bookplot(n,tag,'m(l+,l-)',mll,wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1
                  
      else

c--- lepton+jets-specific plots
      if (plabel(4) .eq. 'ea') then
        ptl=pt(4,p)
        yl=yrap(4,p)
      else
        ptl=pt(7,p)
        yl=yrap(7,p)
      endif
      
c--- pt(lepton)
      call bookplot(n,tag,'pt(lepton)',ptl,wt,wt2,
     & 0d0,200d0,10d0,'log')
      n=n+1

c--- eta(lepton)
      call bookplot(n,tag,'eta(lepton)',yl,wt,wt2,
     & -2.2d0,2.2d0,0.2d0,'lin')
      n=n+1

      endif

c--- these plots are common to dilepton and lepton+jets channels
      call topreconstruct(p,failed,
     & ptt,yt,pttb,ytb,ptttb,yttb,mttb,mwp,mwm)

      if ((failed) .and. (first .eqv. .false.)) then
c--- make sure to increment by the number of plots in the 'else' section
      n=n+15
      
      else
      
c--- pt(t)
      call bookplot(n,tag,'pt(t)',ptt,wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1
      
c--- pt(tbar)
      call bookplot(n,tag,'pt(tbar)',pttb,wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1

c--- pt(t and tbar)
      call bookplot(n,tag,'pt(t and tbar)',ptt,wt,wt2,
     & 0d0,400d0,10d0,'log')
      call bookplot(n,tag,'pt(t and tbar)',pttb,wt,wt2,
     & 0d0,400d0,10d0,'log')
      n=n+1

c--- y(t)
      call bookplot(n,tag,'y(t)',yt,wt,wt2,
     & -2.6d0,2.6d0,0.2d0,'lin')
      n=n+1

c--- y(tbar)
      call bookplot(n,tag,'y(tbar)',ytb,wt,wt2,
     & -2.6d0,2.6d0,0.2d0,'lin')
      n=n+1

c--- y(t and tbar)
      call bookplot(n,tag,'y(t and tbar)',yt,wt,wt2,
     & -2.6d0,2.6d0,0.2d0,'lin')
      call bookplot(n,tag,'y(t and tbar)',ytb,wt,wt2,
     & -2.6d0,2.6d0,0.2d0,'lin')
      n=n+1

c--- pt(t,tbar)
      call bookplot(n,tag,'pt(t,tbar)',ptttb,wt,wt2,
     & -0.01d0,399.99d0,10d0,'log')
      n=n+1

c--- y(t,tbar)
      call bookplot(n,tag,'y(t,tbar)',yttb,wt,wt2,
     & -2.6d0,2.6d0,0.2d0,'lin')
      n=n+1

c--- m(t,tbar)
      call bookplot(n,tag,'m(t,tbar)',mttb,wt,wt2,
     & 0d0,1400d0,25d0,'log')
      n=n+1

c--- m(W+ candidate)
      call bookplot(n,tag,'m(W+ candidate)',mwp,wt,wt2,
     & 0d0,200d0,5d0,'log')
      n=n+1

c--- m(W- candidate)
      call bookplot(n,tag,'m(W- candidate)',mwm,wt,wt2,
     & 0d0,200d0,5d0,'log')
      n=n+1
      endif
      
c--- This section tailored for comparison with CMS PAS TOP-11-013

c--- plots for comparison with arXiv:0907.3090
      if (dilepton) then

c--- pt(l+)
      call bookplot(n,tag,'pt(l+)',pt(4,p),wt,wt2,
     & 0d0,1000d0,50d0,'log')
      n=n+1
c--- pt(l-)
      call bookplot(n,tag,'pt(l-)',pt(7,p),wt,wt2,
     & 0d0,1000d0,50d0,'log')
      n=n+1

c--- eta(l+)
      call bookplot(n,tag,'eta(l+)',yrap(4,p),wt,wt2,
     & -4d0,4d0,0.2d0,'lin')
      n=n+1
c--- eta(l-)
      call bookplot(n,tag,'eta(l-)',yrap(7,p),wt,wt2,
     & -4d0,4d0,0.2d0,'lin')
      n=n+1
      endif

c--- m(l+,b)
      mlb=dsqrt((p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
     &         -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2)
      call bookplot(n,tag,'m(l+,b)',mlb,wt,wt2,
     & 0d0,200d0,5d0,'log')
      n=n+1
c--- m(l-,bb)
      mlb=dsqrt((p(6,4)+p(7,4))**2-(p(6,1)+p(7,1))**2
     &         -(p(6,2)+p(7,2))**2-(p(6,3)+p(7,3))**2)
      call bookplot(n,tag,'m(l-,bb)',mlb,wt,wt2,
     & 0d0,200d0,5d0,'log')
      n=n+1
c--- plots for comparison with arXiv:0907.3090
      

c--- plots for computing ttbar asymmetry
  
      if ((failed) .and. (first .eqv. .false.)) then
c--- make sure to increment by the number of plots in the 'else' section
      n=n+3
      
      else
            
c--- inclusive case
      call bookplot(n,tag,'y(t)-y(tbar)',yt-ytb,wt,wt2,
     & -4d0,4d0,4d0,'lin')
      n=n+1

c--- for bins of |y(t)-y(tbar)| of (0,0.5), (0.5,1), (1,1.5)
      call bookplot(n,tag,'y(t)-y(tbar)',yt-ytb,wt,wt2,
     & -1.5d0,1.5d0,0.5d0,'lin')
      n=n+1

c--- for bins of |y(t)-y(tbar)| of (1.5,4.5)
      call bookplot(n,tag,'y(t)-y(tbar)',yt-ytb,wt,wt2,
     & -4.5d0,4.5d0,3d0,'lin')
      n=n+1

      endif
      
c-- lepton asymmetry
      if (dilepton .eqv. .false.) then
        if (plabel(4) .eq. 'ea') then
        call bookplot(n,tag,'Q(lep)*y(lep)',yrap(4,p),wt,wt2,
     &   -4d0,4d0,4d0,'lin')
        n=n+1
        endif
        if (plabel(7) .eq. 'el') then
        call bookplot(n,tag,'Q(lep)*y(lep)',-yrap(7,p),wt,wt2,
     &   -4d0,4d0,4d0,'lin')
        n=n+1
        endif
      endif
      
c--- plots for computing ttbar asymmetry

c--- invariant mass of two leading jets       
      call bookplot(n,tag,'m(j1,j2)',mj1j2,wt,wt2, 0d0,300d0,5d0,'lin')
      n=n+1
      
c--- compute energy of particle 4 in rest frame of top
      do j=1,4
      plep(j)=p(4,j)
      ptop(j)=p(3,j)+p(4,j)+p(5,j)
      ptoprest(j)=0d0
      enddo
      ptoprest(4)=mt
      call boostx(plep,ptop,ptoprest,pleprest)      
      call bookplot(n,tag,'2*E(4)/mt (top rest frame)',
     & 2d0*pleprest(4)/mt,wt,wt2,0d0,1d0,0.01d0,'lin')
      n=n+1

c--- compute energy of particle 7 in rest frame of top
      do j=1,4
      plep(j)=p(7,j)
      ptop(j)=p(6,j)+p(7,j)+p(8,j)
      ptoprest(j)=0d0
      enddo
      ptoprest(4)=mt
      call boostx(plep,ptop,ptoprest,pleprest)      
      call bookplot(n,tag,'2*E(7)/mt (top rest frame)',
     & 2d0*pleprest(4)/mt,wt,wt2,0d0,1d0,0.01d0,'lin')
      n=n+1
       
c--- single-particle plots
      do j=3,8
        call genplot1(p,j,tag,wt,wt2,n)
      enddo
c--- two-particle plots
      call genplot2(p,3,4,tag,wt,wt2,n)
      call genplot2(p,7,8,tag,wt,wt2,n)
      call genplot2(p,5,6,tag,wt,wt2,n)
c--- three-particle plots
      call genplot3(p,3,4,5,tag,wt,wt2,n)
      call genplot3(p,6,7,8,tag,wt,wt2,n)

c--- additional plots that may be present at NLO
      if (abs(p(9,4)) .gt. 1d-8) then
        call genplot1(p,9,tag,wt,wt2,n)
        call genplot3(p,7,8,9,tag,wt,wt2,n)
      else
        n=n+2
      endif

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return
      end
