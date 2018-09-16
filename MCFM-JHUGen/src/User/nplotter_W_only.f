      subroutine nplotter_W_only(p,wt,wt2,switch)
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
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'jetlabel.f'
      include 'outputflags.f'
      include 'nproc.f'
      real(dp):: p(mxpart,4),wt,wt2,yrap,pt,r,yraptwo,etaraptwo,
     & y3,y4,y5,pt3,pt4,pt5,Re5,y34,eta34
      integer:: switch,n,nplotmax
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,y3,y4,y5,pt5,Re5)
      
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        y3=1d3
        y4=1d3
c--- If there is no NLO jet, these initial y5, pt5 will not pass the cut
        y5=1d3
        pt5=1d3
c--- If Re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        Re5=1d3
        jets=1
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

c--- W rapidity and pseudorapidity
      y34=yraptwo(3,4,p)
      eta34=etaraptwo(3,4,p)

c--- If nproc=1, plot e^+(4). If nproc=6, plot e^-(3).
      if(nproc == 1) then
         y4=yrap(4,p)
         pt4=pt(4,p)
      else
         y3=yrap(3,p)
         pt3=pt(3,p)
      endif
c---      eventpart=4+jets
c---      print*, nproc
      if(jets > 0) then
         pt5=pt(5,p)
         y5=yrap(5,p)
         if(nproc == 1) then
            Re5=R(p,4,5)
         else
            Re5=R(p,3,5)
         endif
      else
        pt5=-1._dp
        y5=1d3
        Re5=1d3
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

       call bookplot(n,tag,'W rapidity',y34,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
       call ebookplot(n,tag,y34,wt)
       n=n+1
       call bookplot(n,tag,'W ps-rap',eta34,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
       n=n+1
      if(nproc == 1) then
         call bookplot(n,tag,'y(lep)',y4,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
         n=n+1
         call bookplot(n,tag,'pt(lep)',pt4,wt,wt2,0._dp,100._dp,2._dp,'lin')
         n=n+1
      else            
         call bookplot(n,tag,'y(lep)',y3,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
         n=n+1
         call bookplot(n,tag,'pt(lep)',pt3,wt,wt2,0._dp,100._dp,2._dp,'lin')
         n=n+1
      endif
      call bookplot(n,tag,'DeltaRe5',Re5,wt,wt2,0._dp,5._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-3._dp,3._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,80._dp,2._dp,'lin')
      n=n+1

  
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
      
