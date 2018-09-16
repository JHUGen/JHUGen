      subroutine nplotter_Ztj(p,wt,wt2,switch)
      implicit none
      include 'types.f'
      
c -- R. Rontsch, 2013-02-28
c -- Taylored for generic observables in pp -> Z(-> e-(p3) e+(p4)) t(p5) d(p6) g(p7) [top not decaying]
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
      include 'plabel.f'
      include 'outputflags.f'
      real(dp):: p(mxpart,4),wt,wt2,ptjet(mxpart),yjet(mxpart),
     & pt,pttwo,yrap,yraptwo,yrapthree,etarap,etaraptwo,etarapthree,
     & tiny
      real(dp):: ptl,pta,ptt,ptj1,ptj2,ptla,
     & yl,ya,yt,yj1,yj2,yla,ytj,
     & ylat,ylaj,etat,etala,etatj,etalat,etalaj,mla,mtj,mlat,mlaj,
     & Dytj,Dylat,Dylaj,Detalat,Detalaj
      real(dp):: twomass,threemass
      integer:: switch,n,nplotmax,j,iorder(mxpart),ijet,hj,sj
      integer tag
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
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        goto 99
      else
c---  Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

c -- start with jet definitions
      ijet=0
       do j=2,7
         if ((plabel(j) == 'pp') .or. (plabel(j) == 'bq')
     &   .or.(plabel(j) == 'ba')) then
           if (p(j,4) > tiny) then
             ijet=ijet+1
             ptjet(ijet)=pt(j,p)
             yjet(ijet)=yrap(j,p)
           endif
         endif
       enddo

c -- arrange jets by hardness      
       call arraysort(ijet,ptjet,iorder)
c -- hj labels hardest jet, sj the next hardest jet      
c       hj=jetindex(iorder(1))
c       sj=jetindex(iorder(2))
       hj=iorder(1)
       sj=iorder(2)
c -- get pt and rapidity of hardest and next-hardest jets
       ptj1=ptjet(hj)
       yj1=yjet(hj)
       if (ijet == 2) then
          ptj2=ptjet(sj)
          yj2=yjet(sj)
       else
          ptj2=0._dp
          yj2=0._dp
       endif
       
       if (ptj2 >= ptj1) then
          write(*,*) ptj1,ptj2
          stop
       endif

c -- pts
       ptl=pt(3,p)
       pta=pt(4,p)
       ptt=pt(5,p)
       ptla=pttwo(3,4,p)

c -- rapidities       
       yl=yrap(3,p)
       ya=yrap(4,p)      
       yt=yrap(5,p)
       yla=yraptwo(3,4,p)
c  -- use hardest jet when looking at multiparticle observables
       ytj=yraptwo(5,hj,p)
       ylat=yrapthree(3,4,5,p)
       ylaj=yrapthree(3,4,hj,p)

c -- pseudorapidities
       etat=etarap(5,p)
       etala=etaraptwo(3,4,p)
       etatj=etaraptwo(5,hj,p)
       etalat=etarapthree(3,4,5,p)
       etalaj=etarapthree(3,4,hj,p)

c -- inv masses
       mla=twomass(3,4,p)
       mtj=twomass(5,hj,p)
       mlat=threemass(3,4,5,p)
       mlaj=threemass(3,4,hj,p)

c --  differences in rapidity/pseudorapidity
       Dytj=yt-yj1
       Dylat=yla-yt
       Dylaj=yla-yj1
       Detalat=etala-etat
       Detalaj=etala-yj1

***************** END DEFINITIONS **************

       


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

c -- One-particle plots
c -- lepton plots
      call bookplot(n,tag,'pt(e-)',ptl,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(e-)',yl,wt,wt2,
     & -2.2_dp,2.2_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(e+)',pta,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(e+)',ya,wt,wt2,
     & -2.2_dp,2.2_dp,0.2_dp,'lin')
      n=n+1     

c -- top plots
      call bookplot(n,tag,'pt(t)',ptt,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'eta(t)',etat,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'y(t)',yt,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1     

c -- jet plots
      call bookplot(n,tag,'pt(j1)',ptj1,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'y(j1)',yj1,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'pt(j2)',ptj2,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'y(j2)',yj2,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1     

c -- Two particle plots
      call bookplot(n,tag,'pt(la)',ptla,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'m(Z)',mla,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1     
      call bookplot(n,tag,'m(Z)',mla,wt,wt2,0._dp,200._dp,5._dp,'lin')
      n=n+1     
      call bookplot(n,tag,'m(tj)',mtj,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'y(Z)',yla,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'eta(Z)',etala,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'y(tj)',ytj,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'eta(tj)',etatj,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'Delta y(tj)',Dytj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1     
      
c -- Three particle plots
      call bookplot(n,tag,'m(lat)',mlat,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'m(laj)',mlaj,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1     
      call bookplot(n,tag,'y(lat)',ylat,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'eta(lat)',etalat,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'y(laj)',ylaj,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'eta(laj)',etalaj,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'Delta y(la,t)',Dylat,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'Delta eta(la,t)',Detalat,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'Delta y(la,j)',Dylaj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1     
      call bookplot(n,tag,'Delta eta(la,j)',Detalaj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
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
