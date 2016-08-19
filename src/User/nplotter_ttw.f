      subroutine nplotter_ttw(p,wt,wt2,switch)
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
      include 'plabel.f'
      include 'npart.f'
      include 'outputflags.f'
      double precision p(mxpart,4),wt,wt2,pt,tiny,HTjet,MET,
     & etmiss,vecmet(4),binindex,METHTdoublebin
      integer switch,n,nplotmax,j
      character*4 tag
      parameter(tiny=1d-8)
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

c--- missing ET
      MET=etmiss(p,vecmet)
c--- HT(jet) - scalar sum of jet pt
      HTjet=0d0
      do j=3,npart+2-switch
        if ((plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'bq')
     & .or. (plabel(j) .eq. 'ba')) then
          HTjet=HTjet+pt(j,p)
        endif
      enddo
       
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

      call bookplot(n,tag,'HTjet',HTjet,wt,wt2,0d0,800d0,20d0,'lin')
      n=n+1
      call bookplot(n,tag,'MET',MET,wt,wt2,0d0,200d0,10d0,'lin')
      n=n+1
      
c--- This section tailored for comparison with CMS PAS SUS-11-020   

c--- only fill histograms if both b-quark and anti-b-quark are present
      if ((first) .or. ((p(5,4) .gt. tiny).and.(p(6,4) .gt. tiny))) then
    
      call bookplot(n,tag,'HTjet with 2b',HTjet,wt,wt2,
     & 0d0,800d0,20d0,'lin')
      n=n+1
      call bookplot(n,tag,'MET with 2b',MET,wt,wt2,0d0,200d0,10d0,'lin')
      n=n+1
      
      if ((first) .or. ((MET .gt. 30d0) .and. (HTjet .gt. 80d0))) then
      call bookplot(n,tag,'SR1',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 120d0) .and. (HTjet .gt. 200d0))) then
      call bookplot(n,tag,'SR3',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      call ebookplot(n,tag,0.5d0,wt)
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 50d0) .and. (HTjet .gt. 200d0))) then
      call bookplot(n,tag,'SR4',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      call ebookplot(n,tag,0.5d0,wt)
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 50d0) .and. (HTjet .gt. 320d0))) then
      call bookplot(n,tag,'SR5',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      call ebookplot(n,tag,0.5d0,wt)
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 120d0) .and. (HTjet .gt. 320d0))) then
      call bookplot(n,tag,'SR6',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      call ebookplot(n,tag,0.5d0,wt)
      endif
      n=n+1

c--- double-binned histogram in MET and HT
      binindex=METHTdoublebin(MET,HTjet)
      call bookplot(n,tag,'MET/HT double bin',binindex,wt,wt2,
     & -0.5d0,99.5d0,1d0,'lin')
      n=n+1     
           
      else

c--- increment n by the # of histograms skipped      
      n=n+8
      
      endif
         
c--- This section tailored for comparison with CMS PAS SUS-11-020       

c--- This section tailored for comparison with CMS PAS SUS-11-010   
      
      if ((first) .or. ((MET .gt. 120d0) .and. (HTjet .gt. 400d0))) then
      call bookplot(n,tag,'Region 1',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 50d0) .and. (HTjet .gt. 400d0))) then
      call bookplot(n,tag,'Region 2',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 120d0) .and. (HTjet .gt. 200d0))) then
      call bookplot(n,tag,'Region 3',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
              
      if ((first) .or. ((MET .gt. 100d0) .and. (HTjet .gt. 80d0))) then
      call bookplot(n,tag,'Region 4',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      endif
      n=n+1
              
         
c--- This section tailored for comparison with CMS PAS SUS-11-010       

              
c--- single-particle plots
      do j=3,10
        call genplot1(p,j,tag,wt,wt2,n)
      enddo
c--- two-particle plots
      call genplot2(p,3,4,tag,wt,wt2,n)
      call genplot2(p,7,8,tag,wt,wt2,n)
      call genplot2(p,9,10,tag,wt,wt2,n)
      call genplot2(p,5,6,tag,wt,wt2,n)
c--- three-particle plots
      call genplot3(p,3,4,5,tag,wt,wt2,n)
      call genplot3(p,6,7,8,tag,wt,wt2,n)

c--- additional plots that may be present at NLO       
      if (abs(p(11,4)) .gt. 1d-8) then
        call genplot1(p,11,tag,wt,wt2,n)
      else
        n=n+1
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
