      subroutine nplotter_auto(p,wt,wt2)
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
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'outputflags.f'
      include 'mcfmplotinfo.f'
      double precision p(mxpart,4),wt,wt2,tiny
      integer n,nplotmax,j,m,i1,i2,i3,i4,ilomomenta
      character*4 tag
      common/ilomomenta/ilomomenta
      common/nplotmax/nplotmax
      parameter(tiny=1d-8)
      logical, save::first=.true.
ccccc!$omp threadprivate(first,/nplotmax/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c---- Initialize histograms and determine the number of momenta
c---- present at leading order by inspecting plabel
        tag='book'
      else
c--- Add event in histograms
        tag='plot'
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************
       
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

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

c--- single-particle plots
      do j=3,ilomomenta
        call autoplot1(p,j,tag,wt,wt2,n)
      enddo
c--- two-particle plots
      
      j=1
      do while (mcfmplotinfo(j) > 0)
        m=mcfmplotinfo(j)
        if     (m .lt. 10) then
c--- one-particle plots
          i1=m
          if (i1 .eq. 0) i1=10   ! special code: 0 -> 10
          call autoplot1(p,i1,tag,wt,wt2,n)
        elseif (m .lt. 100) then
c--- two-particle plots
          i1=m/10
          i2=mod(m,10)
          if (i2 .eq. 0) i2=10   ! special code: 0 -> 10
          call autoplot2(p,m,i1,i2,tag,wt,wt2,n)
        elseif (m .lt. 1000) then
c--- three-particle plots
          i1=m/100
          i2=(m-i1*100)/10
          i3=mod(m,10)
          if (i3 .eq. 0) i3=10   ! special code: 0 -> 10
          call autoplot3(p,m,i1,i2,i3,tag,wt,wt2,n)
        elseif (m .lt. 10000) then
c--- four-particle plots
          i1=m/1000
          i2=(m-i1*1000)/100
          i3=(m-i1*1000-i2*100)/10
          i4=mod(m,10)
          if (i4 .eq. 0) i4=10   ! special code: 0 -> 10
          call autoplot4(p,m,i1,i2,i3,i4,tag,wt,wt2,n)
        else
          write(6,*) 'Unforeseen plot in nplotter_auto.f'
          stop
        endif
        j=j+1
      enddo
         
c--- additional plots that may be present at NLO       
      if ((abs(p(ilomomenta+1,4)).gt.tiny) .or. (tag.eq.'book')) then
        call autoplot1(p,ilomomenta+1,tag,wt,wt2,n)
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
