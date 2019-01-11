      subroutine nplotter_gamgam(p,wt,wt2,switch,nd)
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
!-- nd determines whether we are analysing a photon dipole and hence have
!---to rescale accordingly
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      double precision p(mxpart,4),wt,wt2
      double precision yrap,pt,r,yraptwo
      double precision y3,y4,y34,y5,pt3,pt4,pt5,r45,r35,r34,s34,m34
      integer switch,n,nplotmax
      character*4 tag
      integer nd
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/)

  
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
        pt3=1d3
        pt4=1d3
        y34=1d3
c---Intiailise jet 
        y5=1d3
        pt5=1d3
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        r34=1d3
        r35=1d3
        r45=1d3
        jets=nqcdjets
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
c--- Photons order based on pt
      if(dabs(pt(3,p)).gt.dabs(pt(4,p))) then 
         pt3=pt(3,p)
         pt4=pt(4,p)
         y3=0.5d0*(yrap(3,p)+yrap(4,p))
         if(jets.gt.0) then 
             r35=R(p,3,5)
             r45=R(p,4,5)
          endif
       else
          pt3=pt(4,p)
          pt4=pt(3,p)
          y3=0.5d0*(yrap(3,p)+yrap(4,p))
          if(jets.gt.0) then 
             r35=R(p,4,5)
             r45=R(p,3,5)
          endif
       endif

      !y34
      y34=yraptwo(3,4,p)
      

      ! Jet
      if(jets.gt.0) then 
         pt5=pt(5,p)
         y5=yrap(5,p)
      else ! put out of range
         pt5=1d7
         y5=1d7
         r45=1d7
         r35=1d7
      endif
      
      ! R_ij
      r34=R(p,3,4)

      ! m_gamma,gamma
      s34=2d0*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
  
      m34=dsqrt(s34)

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
      call bookplot(n,tag,'1/2(y3+y4)',y3,wt,wt2,
     & -6d0,6d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0d0,200d0,5d0,'log')
      n=n+1

      call bookplot(n,tag,'pt3_2',pt3,wt,wt2,0d0,80d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'pt4_2',pt4,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
                
      call bookplot(n,tag,'y34',y34,wt,wt2,-6d0,6d0,0.5d0,'lin')
      n=n+1
     
      call bookplot(n,tag,'DeltaR34',r34,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
     
!    m34 
      call bookplot(n,tag,'m34',m34,wt,wt2,100d0,150d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,0d0,300d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'y5',y5,wt,wt2,-3d0,3d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0d0,200d0,2d0,'log')
      n=n+1
   
      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR45',r45,wt,wt2,0d0,5d0,0.1d0,'lin')
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
      
      

