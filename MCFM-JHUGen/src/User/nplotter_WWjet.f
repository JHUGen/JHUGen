      subroutine nplotter_WW_jet(p,wt,wt2,switch)
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
      include 'outputflags.f'
      double precision p(mxpart,4),wt,wt2,pt,pord(mxpart,4),
     & pt7,pt8,mll,delphi,etvec(4),Etll,MET,mtrans
      integer switch,n,nplotmax
      character*4 tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,mll,delphi,MET,mtrans)


************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        mll=1d7
        delphi=1d7
        MET=1d7
        mtrans=1d7
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

c--- dilepton invariant mass
      mll=
     & +(p(4,4)+p(5,4))**2
     & -(p(4,1)+p(5,1))**2
     & -(p(4,2)+p(5,2))**2
     & -(p(4,3)+p(5,3))**2
      mll=sqrt(max(0d0,mll))
     
c--- dilepton azimuthal separation
      delphi=(p(4,1)*p(5,1)+p(4,2)*p(5,2))
     &       /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (delphi .gt. +0.9999999D0) delphi=+1d0
      if (delphi .lt. -0.9999999D0) delphi=-1d0
      delphi=dacos(delphi)

c--- missing ET
      etvec(:)=p(3,:)+p(6,:)
      MET=sqrt(max(0d0,etvec(1)**2+etvec(2)**2))

c--- transverse mass
      Etll=sqrt(max(0d0,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
      mtrans=MET+Etll
     
      pord(:,:)=p(:,:)
c--- order jets by pt if there are two present
      if (abs(p(8,4)) .gt. 1d-8) then
        pt7=pt(7,p)
        pt8=pt(8,p)
        if (pt8 .gt. pt7) then
          pord(7,:)=p(8,:)
          pord(8,:)=p(7,:)
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

c--- usual plots for particle 3 (neutrino)
      call autoplot1(p,3,tag,wt,wt2,n)

c--- usual plots for particle 4 (e+)
      call autoplot1(p,4,tag,wt,wt2,n)

c--- usual plots for particle 5 (e-)
      call autoplot1(p,5,tag,wt,wt2,n)

c--- usual plots for particle 6 (antineutrino)
      call autoplot1(p,6,tag,wt,wt2,n)

c--- usual plots for particle 7 (highest pt jet)
      call autoplot1(p,7,tag,wt,wt2,n)

c--- usual plots for W+
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for W-
      call autoplot2(p,56,5,6,tag,wt,wt2,n)

c--- dilepton invariant mass
      call bookplot(n,tag,'mll',mll,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
      
c--- dilepton azimuthal separation
      call bookplot(n,tag,'delphi',delphi,wt,wt2,0d0,3.14d0,0.1d0,'lin')
      n=n+1

c--- missing ET
      call bookplot(n,tag,'MET',MET,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
      
c--- transverse mass
      call bookplot(n,tag,'mtrans',mtrans,wt,wt2,0d0,1000d0,20d0,'log')
      n=n+1
      
c--- additional plots that may be present at NLO
      if (abs(p(8,4)) .gt. 1d-8) then
        call autoplot1(pord,8,tag,wt,wt2,n)
        call autoplot2(pord,78,7,8,tag,wt,wt2,n)
      else
        n=n+5
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
