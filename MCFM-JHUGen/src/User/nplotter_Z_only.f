      subroutine nplotter_Z_only(p,wt,wt2,switch)
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
      include 'jetlabel.f'
      include 'outputflags.f'
      double precision p(mxpart,4),wt,wt2
      double precision yrap,pt,yraptwo,pttwo,r
c---  Z->e+e-(31) or b bbar(33): both measured, rapidities and momenta of 3 and 4 can
c---  be calculated, also the invariant mass m34
      double precision y3,y4,y5,y34,pt3,pt4,pt5,pt34,m34,r35
      double precision costheta,p3(4),p4(4),p34(4)
      integer switch,n,nplotmax
      character*4 tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,y4,y5,y34,pt3,pt4,pt5,pt34,m34,r35)


************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        y3=1d3
        y4=1d3
        y5=1d3
        y34=1d3
        pt3=0d0
        pt4=0d0
        pt5=1d3
        pt34=0d0
        m34=0d0
        r35=1d3
        jets=1
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

      y3=yrap(3,p)
      y4=yrap(4,p)
      y34=yraptwo(3,4,p)
      pt3=pt(3,p)
      pt4=pt(4,p)
      pt34=pttwo(3,4,p)
      m34=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)

      if(jets .gt. 0) then
         pt5=pt(5,p)
         y5=yrap(5,p)
         r35=R(p,3,5)
      else
         pt5=-1d0
         y5=1d3
         r35=1d3
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
           
      call bookplot(n,tag,'y3',y3,wt,wt2,-6d0,6d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'y4',y4,wt,wt2,-6d0,6d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-6d0,6d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0d0,80d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,0d0,80d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt34',pt34,wt,wt2,0d0,50d0,2d0,'lin')
      n=n+1
      call bookplot(n, tag,'m34',m34,wt,wt2,70d0,110d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-3.2d0,3.2d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0d0,100d0,2d0,'lin')
      n=n+1
      
c--- compute lepton asymmetry as a function of m34  
c--- (see for example Eq.(3) of PLB718 (2013) 752)
      p3(:)=p(3,:)
      p4(:)=p(4,:)
      p34(:)=p(3,:)+p(4,:)
      costheta=p34(3)/abs(p34(3))
     & *((p3(4)+p3(3))*(p4(4)-p4(3))-(p3(4)-p3(3))*(p4(4)+p4(3)))
     & /m34/sqrt(m34**2+p34(1)**2+p34(2)**2)
c--- these histograms must be kept
      if ((costheta .gt. 0d0) .or. (tag .eq. 'book')) then
        call bookplot(n, tag,'m34 forward lepton',
     &   m34,wt,wt2,40d0,200d0,5d0,'lin')
      endif
      n=n+1
      if ((costheta .le. 0d0) .or. (tag .eq. 'book')) then
        call bookplot(n, tag,'m34 backward lepton',
     &   m34,wt,wt2,40d0,200d0,5d0,'lin')
      endif
      n=n+1
c--- now compute asymmetry - histograms n+1 and n+2 are only temporary
      if (tag .eq. 'book') then
        call bookplot(n, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40d0,200d0,5d0,'lin')
        call bookplot(n+1, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40d0,200d0,5d0,'lin')
        call bookplot(n+2, tag,'lepton FB asymmetry',
     &   m34,wt,wt2,40d0,200d0,5d0,'lin')
      endif
      call mopera(n-2,'-',n-1,n,1d0,1d0)
      call mopera(n-2,'+',n-1,n+1,1d0,1d0)
c--- this is the histogram (n) we will keep
      call mopera(n,'/',n+1,n,1d0,1d0)
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
      
