      subroutine nplotter_VHWW(p,wt,wt2,switch,nd)
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
      real(dp):: p(mxpart,4),wt,wt2,m3456,pt34,pttwo
      integer:: switch,n,nplotmax
      real(dp) :: mbb,mVbb,ptbb,ptbh,ptbs,ptV
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
      real(dp) :: pt5,pt6,pt56,pt,pt3,pt4,mtW,mVHtrans,mymet
      real(dp) :: pl1,pl2,pl3,pl4
      integer nj,nbj,noj
      common/observables_VHWW/mVHtrans,mymet,pl1,pl2,pl3,pl4
!$omp threadprivate(/observables_VHWW/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

!      write(6,*) nj,nd

      nj=jets

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


!=====this is MVH, defined if two b-jets are presente
 !     if(nbj>=2) then 
 !        m3456=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
 !    &        -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
 !    &        -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
 !    &        -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
 !        m3456=sqrt(max(m3456,zip))
 !     else
 !        m3456=-1._dp
 !     endif

      pt3=pt(3,p)
      pt4=pt(4,p)
      pt34=pttwo(3,4,p)
!      pt5=pt(5,p)
!      pt6=pt(6,p) 
!      pt56=pttwo(5,6,p)

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

      call bookplot(n,tag,'njets ',
     &     real(nj,dp),wt,wt2,-0.5_dp,5.5_dp,1._dp,'lin')
      n=n+1
  
      call bookplot(n,tag,'mVHtrans',
     & mVHtrans,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1

      
      call bookplot(n,tag,'missing ET',
     & mymet,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt lep 1',
     &     pl1,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt lep 2',
     &     pl2,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt lep 3',
     &     pl3,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt lep 4',
     &     pl4,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      
c---  usual plots for 3+4
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for 5+6
      call autoplot2(p,56,5,6,tag,wt,wt2,n)
      call autoplot2(p,78,7,8,tag,wt,wt2,n)

c--- usual plots for 5+6+7+8
      call autoplot4(p,5678,5,6,7,8,tag,wt,wt2,n)
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
