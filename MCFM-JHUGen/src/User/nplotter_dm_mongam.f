      subroutine nplotter_dm_mongam(p,wt,wt2,switch,nd)
      implicit none
      include 'types.f'
!==== mono-jet Nplotter
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
      integer:: i,nd
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: etarap,pt,etaj1,ptj1,pt34,pttwo,m34,m345
      real(dp):: ptga,etaga
      integer:: switch,n,nplotmax
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/)
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************


c--- Initialize dummy values for all quantities that could be plotted
      ptj1=-1._dp
      ptga=-1._dp 
      etaga=99._dp
      etaj1=99._dp
      m34=-1._dp 
      m345=-1._dp 
      
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


!======= Missing ET 
      pt34=pttwo(3,4,p)
!====== DM invariant mass (NOT directly observable, but useful for theory) 
      m34=0._dp 
!======= DM + jet invariant mass (NOT directly observable, but useful for theory) 
      m345=0._dp 
      do i=1,3
         m34=m34-(p(3,i)+p(4,i))**2 
         m345=m345-(p(3,i)+p(4,i)+p(5,i))**2
      enddo
      m34=m34+(p(3,4)+p(4,4))**2 
      m345=m345+(p(3,4)+p(4,4)+p(5,4))**2 
      m345=sqrt(m345)
      m34=sqrt(m34) 
      

!========== phpton pt 
      ptga=pt(5,p) 
      etaga=etarap(5,p) 
 

!====== jet pt and eta
      if (jets == 1) then 
         ptj1=pt(6,p)
         etaj1=etarap(6,p) 
      else  !=== out of plotting range
         ptj1=-1._dp 
          etaj1=99._dp 
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

      call bookplot(n,tag,'Missing ET',pt34,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'DM inv mass',m34,wt,wt2,
     &              0._dp,1200._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Photon pt',ptga,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Photon pt lin',ptga,wt,wt2,
     &              0._dp,2000._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'Photon eta',etaga,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt lin',ptj1,wt,wt2,
     &              0._dp,2000._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 eta',etaj1,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
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
      
