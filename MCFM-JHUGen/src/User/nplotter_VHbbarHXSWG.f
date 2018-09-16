      subroutine nplotter_VHbbarHXSWG(p,wt,wt2,switch,nd)
      implicit none
      include 'types.f'
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
c--- switch:  an integer:: equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      include 'interference.f'
      include 'jetlabel.f'
      integer nd
      real(dp):: p(mxpart,4),wt,wt2,m3456,pttwo,yrap,yraptwo
      integer:: switch,n,nplotmax,iloop
      real(dp) :: mbb,mVbb,ptbb,ptbh,ptbs,ptV,ptH,yH,y3,y4,pt3,pt4
      integer tag
      character*9 range
      logical failedbcuts
      logical, save::first=.true.
      common/nplotmax/nplotmax
      real(dp) :: pt5,pt6,pt56,pt,ayrap,mtW
      integer nj,nbj,noj,nbq,nba
      common/njetsVH/nj,nbj,noj
      common/observables_VHbb/ptV,mbb,mVbb,ptbb,ptbh,ptbs,mtW
!$omp threadprivate(/njetsVH/,/observables_VHbb/) 

ccccc!$omp threadprivate(first,/nplotmax/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

      pt3=pt(3,p)
      pt4=pt(4,p)
      ptV=pttwo(3,4,p)

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
        call bookfill(tag,p,wt/real(itmx,dp))  
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
 
      failedbcuts=.false.
      call getbs(p,nbq,nba)
      if ((nbq <= 0) .or. (nba <=0)) then
        failedbcuts=.true.
      endif
      
      if (pt(nbq,p) < 25._dp) failedbcuts=.true.
      if (pt(nba,p) < 25._dp) failedbcuts=.true.
      if (ayrap(nbq,p) > 2.5_dp) failedbcuts=.true.
      if (ayrap(nba,p) > 2.5_dp) failedbcuts=.true.
      
      if (failedbcuts) then
c--- put variables out of range of plots
        ptV=-1d0
        ptH=-1d0
        yH=99d0
        pt3=-1d0
        pt4=-1d0
        y3=99d0
        y4=99d0
      else
        ptH=pttwo(nbq,nba,p)
        yH=yraptwo(nbq,nba,p)
        pt3=pt(3,p)
        pt4=pt(4,p)
        y3=yrap(3,p)
        y4=yrap(4,p)
      endif
 
      call bookplot(n,tag,'ptV',ptV,wt,wt2,0._dp,500._dp,5._dp,'log')
      n=n+1
      
      do iloop=0,3
        if (iloop == 0) range='inclusive'
        if (iloop == 1) range='[0,150\  '
        if (iloop == 2) range='[150,250\'
        if (iloop == 3) range='> 250    '
        if (
     & (iloop == 1 .and. (ptV > 150._dp)) .or.
     & (iloop == 2 .and. ((ptV < 150._dp) .or. (ptV > 250._dp))) .or.
     & (iloop == 3 .and. (ptV < 250._dp)) ) then
          if (first .eqv. .false.) then
            n=n+6
            cycle
          endif
        endif
        call bookplot(n,tag,'ptH, ptV '//range,
     &    ptH,wt,wt2,0._dp,500._dp,5._dp,'log')
        n=n+1
        call bookplot(n,tag,'yH, ptV '//range,
     &    yH,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt3, ptV '//range,
     &    pt3,wt,wt2,0._dp,500._dp,5._dp,'log')
        n=n+1
        call bookplot(n,tag,'y3, ptV '//range,
     &    y3,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt4, ptV '//range,
     &    pt3,wt,wt2,0._dp,500._dp,5._dp,'log')
        n=n+1
        call bookplot(n,tag,'y4, ptV '//range,
     &    y4,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
        n=n+1
      enddo      
      
c--- usual plots for 3+4
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for 3+4+5+6
      call autoplot4(p,3456,3,4,5,6,tag,wt,wt2,n)

!      if(nbj>=2) then 
c--- usual plots for 5+6
!      call autoplot2(p,56,5,6,tag,wt,wt2,n)

c--- usual plots for 3+4+5+6
!      call autoplot4(p,3456,3,4,5,6,tag,wt,wt2,n)
!      endif
c--- additional plots that may be present at NLO       
!      if (abs(p(7,4)) > 1.e-8_dp) then
!        call autoplot1(p,7,tag,wt,wt2,n)
!      else
!        n=n+2
!      endif

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
