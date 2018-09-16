      subroutine nplotter_4ftwdk(p,wt,wt2,switch)
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
      include 'outputflags.f'
      include 'nqcdjets.f'
      real(dp):: p(mxpart,4),wt,wt2,pt
      real(dp):: tiny,swap,costheta,ylight
      integer:: switch,n,nplotmax,j
      integer tag
      logical:: failed
      parameter(tiny=1.e-8_dp)
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms only
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

c--- Jets have already been reordered so that bottom, anti-bottom
c--- quarks are in the correct positions; just have to order jets in
c--- positions 7 and 8 according to pt
      if (p(8,4) > tiny) then
        if (pt(8,p) > pt(7,p)) then
          do j=1,4
          swap=p(7,j)
          p(7,j)=p(8,j)
          p(8,j)=swap
          enddo
        endif
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

c--- fill variables that require top reconstruction
      call singletopreconstruct(p,failed,costheta,ylight)

      if ((failed) .and. (first .eqv. .false.)) then
c--- make sure to increment by the number of plots in the 'else' section
      n=n+2
      
      else
      
c--- cos(theta*)
      call bookplot(n,tag,'cos(theta*)',costheta,wt,wt2,
     & -1._dp,1._dp,0.1_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'|y(light)|',abs(ylight),wt,wt2,
     & 0._dp,5._dp,0.2_dp,'lin')
      n=n+1
      
      endif
      
c--- single-particle plots
      do j=3,7
        if ((abs(p(j,4)) > tiny) .or. (first)) then
          call genplot1(p,j,tag,wt,wt2,n)
        else
          n=n+2
        endif
      enddo
c--- two-particle plots
      call genplot2(p,3,4,tag,wt,wt2,n)
      if (((abs(p(5,4)) > tiny) .and. (abs(p(6,4)) > tiny))
     &    .or. (first)) then
        call genplot2(p,5,6,tag,wt,wt2,n)
      else
        n=n+3
      endif
c--- three-particle plots
      if ((abs(p(5,4)) > tiny) .or. (first)) then
        call genplot3(p,3,4,5,tag,wt,wt2,n)
      else
        n=n+3
      endif

c--- additional plots that may be present at NLO       
      if ((abs(p(8,4)) > tiny) .or. (first)) then
        call genplot1(p,8,tag,wt,wt2,n)
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
      
      



