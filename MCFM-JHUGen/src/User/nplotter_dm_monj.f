      subroutine nplotter_dm_monj(p,wt,wt2,switch)
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
      integer:: i5,i6,i7,nu,i
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: etarap,pt,
     & etaj1,ptj1,etaj2,ptj2,etaj3,mjj,delr,getet,
     & rar,deleta,delfi,pt5,pt6,pt7,tmp5(4),tmp6(4),tmp7(4),oldpt(5:7),
     & sumeta,etastar,pt34,pttwo,m34,m345
      integer:: switch,n,nplotmax
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,m34,m345)
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ptj1=-1._dp
      ptj2=-1._dp
      etaj1=99._dp
      etaj2=99._dp
      etaj3=99._dp
      mjj=-1._dp
      delr=-1._dp
      deleta=-1._dp
      etastar=99._dp
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
      

     
c--- BEGIN: order jets according to pt
      pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
      if (jets > 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
      if (jets > 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))      
      i5=5
      i6=6
      i7=7
      oldpt(5)=pt5
      if (jets > 1) oldpt(6)=pt6
      if (jets > 2) oldpt(7)=pt7
c--- sort for 2 jets 
      if (jets == 2) then          
        if (pt6 > pt5) then
          i5=6
          i6=5
        endif
      endif
c--- sort for 3 jets 
      if (jets == 3) then
        if ((pt5 > pt6) .and. (pt5 > pt7)) then
           i5=5
          if (pt6 > pt7) then
            i6=6
            i7=7
          else
            i6=7
            i7=6
          endif
        endif
        if ((pt6 > pt5) .and. (pt6 > pt7)) then
           i5=6
          if (pt5 > pt7) then
            i6=5
            i7=7
          else
            i6=7
            i7=5
          endif
        endif
        if ((pt7 > pt5) .and. (pt7 > pt6)) then
           i5=7
          if (pt5 > pt6) then
            i6=5
            i7=6
          else
            i6=6
            i7=5
          endif
        endif
      endif
c--- perform exchange
      do nu=1,4
        tmp5(nu)=p(i5,nu)
        tmp6(nu)=p(i6,nu)
        tmp7(nu)=p(i7,nu)
      enddo
      do nu=1,4
        p(5,nu)=tmp5(nu)
        p(6,nu)=tmp6(nu)
        p(7,nu)=tmp7(nu)
      enddo
      pt5=oldpt(i5)
      if (jets > 1) pt6=oldpt(i6)
      if (jets > 2) pt7=oldpt(i7)
c--- END: ordering

c--- Calculate quantities to plot
      etaj1=etarap(5,p)
      ptj1=pt(5,p)
      if (jets >= 2) then
        etaj2=etarap(6,p)
        ptj2=pt(6,p)
        deleta=abs(etaj1-etaj2)
      endif
      if (jets >= 3) then
        etaj3=etarap(7,p)
        if     (abs(etaj1-etaj3) > 
     &      max(abs(etaj1-etaj2),abs(etaj2-etaj3))) then
          deleta=abs(etaj1-etaj3)
          sumeta=etaj1+etaj3
          etaj3=etaj2
        elseif (abs(etaj2-etaj3) > 
     &      max(abs(etaj1-etaj2),abs(etaj1-etaj3))) then
          deleta=abs(etaj2-etaj3)
          sumeta=etaj2+etaj3
          etaj3=etaj1
        else
          deleta=abs(etaj1-etaj2)
          sumeta=etaj1+etaj2
        endif
      endif

      rar=(p(5,1)*p(6,1)+p(5,2)*p(6,2))
     & /sqrt((p(5,1)**2+p(5,2)**2)*(p(6,1)**2+p(6,2)**2))
      if (rar<-1._dp) then 
         delfi=pi
      elseif (rar>1._dp) then 
         delfi=0._dp
      else
         delfi=acos(rar)
      endif
      delr=sqrt(deleta**2+delfi**2)

      if (jets > 1) then
      mjj=sqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &         -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
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
      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt lin',ptj1,wt,wt2,
     &              0._dp,2000._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 eta',etaj1,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      if(jets>1) then 
      call bookplot(n,tag,'Jet 2 eta',etaj2,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt',ptj2,wt,wt2,
     &              0._dp,2000._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt lin',ptj2,wt,wt2,
     &              0._dp,2000._dp,50._dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 eta',etaj2,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'(jet 1,jet 2) invariant mass',mjj,wt,wt2,
     &              0._dp,425._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'delr',delr,wt,wt2,
     &              0.35_dp,4.85_dp,0.15_dp,'lin')
      n=n+1
      call bookplot(n,tag,'deleta',deleta,wt,wt2,
     &              0._dp,4._dp,0.4_dp,'lin')
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
      
