      subroutine nplotter_Wjets(p,wt,wt2,switch)
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
      include 'nproc.f'
      integer i5,i6,i7,nu
      double precision p(mxpart,4),wt,wt2
      double precision etarap,pt,
     & etaj1,ptj1,etaj2,ptj2,etaj3,ptj3,mjj,delr,getet,
     & rar,deleta,delfi,pt5,pt6,pt7,tmp5(4),tmp6(4),tmp7(4),oldpt(5:7),
     & sumeta,etastar,y34,eta34,pt34,y3,pt3,y4,pt4,yrap,yraptwo,
     & etaraptwo,pttwo
      integer switch,n,nplotmax
      character*4 tag
      common/nplotmax/nplotmax
      logical, save::first=.true.
ccccc!$omp threadprivate(first,/nplotmax/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ptj1=-1d0
      ptj2=-1d0
      ptj3=-1d0
      etaj1=99d0
      etaj2=99d0
      etaj3=99d0
      mjj=-1d0
      delr=-1d0
      deleta=-1d0
      etastar=99d0
      
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

c--- W rapidity and pseudorapidity
      y34=yraptwo(3,4,p)
      eta34=etaraptwo(3,4,p)
      pt34=pttwo(3,4,p)

c--- If nproc=11, plot e^+(4). If nproc=16, plot e^-(3).
      if(nproc .eq. 11) then
         y4=yrap(4,p)
         pt4=pt(4,p)
      else
         y3=yrap(3,p)
         pt3=pt(3,p)
      endif

c--- BEGIN: order jets according to pt
      pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
      if (jets .gt. 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
      if (jets .gt. 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))      
      i5=5
      i6=6
      i7=7
      oldpt(5)=pt5
      if (jets .gt. 1) oldpt(6)=pt6
      if (jets .gt. 2) oldpt(7)=pt7
c--- sort for 2 jets 
      if (jets .eq. 2) then          
        if (pt6 .gt. pt5) then
          i5=6
          i6=5
        endif
      endif
c--- sort for 3 jets 
      if (jets .eq. 3) then
        if ((pt5 .gt. pt6) .and. (pt5 .gt. pt7)) then
           i5=5
          if (pt6 .gt. pt7) then
            i6=6
            i7=7
          else
            i6=7
            i7=6
          endif
        endif
        if ((pt6 .gt. pt5) .and. (pt6 .gt. pt7)) then
           i5=6
          if (pt5 .gt. pt7) then
            i6=5
            i7=7
          else
            i6=7
            i7=5
          endif
        endif
        if ((pt7 .gt. pt5) .and. (pt7 .gt. pt6)) then
           i5=7
          if (pt5 .gt. pt6) then
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
      if (jets .gt. 1) pt6=oldpt(i6)
      if (jets .gt. 2) pt7=oldpt(i7)
c--- END: ordering

c--- Calculate quantities to plot
      etaj1=etarap(5,p)
      ptj1=pt(5,p)
      if (jets .ge. 2) then
        etaj2=etarap(6,p)
        ptj2=pt(6,p)
        deleta=abs(etaj1-etaj2)
      endif
      if (jets .ge. 3) then
        etaj3=etarap(7,p)
        ptj3=pt(7,p)
      if     (abs(etaj1-etaj3) .gt. 
     &      max(abs(etaj1-etaj2),abs(etaj2-etaj3))) then
          deleta=abs(etaj1-etaj3)
        sumeta=etaj1+etaj3
        etaj3=etaj2
      elseif (abs(etaj2-etaj3) .gt. 
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
      if (rar.lt.-1d0) then 
         delfi=pi
      elseif (rar.gt.1d0) then 
         delfi=0d0
      else
         delfi=dacos(rar)
      endif
      delr=dsqrt(deleta**2+delfi**2)

      if (jets .gt. 1) then
      mjj=dsqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &         -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
      endif
            
      if ((jets .eq. 3) .and. (deleta .ge. 2d0)) then
        etastar=abs(etaj3-sumeta/2d0)      
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

      call bookplot(n,tag,'W rapidity',y34,wt,wt2,-5d0,5d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'W ps-rap',eta34,wt,wt2,-5d0,5d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'W pt',pt34,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'W pt',pt34,wt,wt2,0d0,1200d0,50d0,'log')
      n=n+1
      if(nproc .eq. 11) then
       call bookplot(n,tag,'y(lep)',y4,wt,wt2,-5d0,5d0,0.4d0,'lin')
       n=n+1
       call bookplot(n,tag,'pt(lep)',pt4,wt,wt2,0d0,200d0,4d0,'lin')
       n=n+1
      else          
       call bookplot(n,tag,'y(lep)',y3,wt,wt2,-5d0,5d0,0.4d0,'lin')
       n=n+1
       call bookplot(n,tag,'pt(lep)',pt3,wt,wt2,0d0,200d0,4d0,'lin')
       n=n+1
      endif

      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0d0,1200d0,50d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt log',ptj1,wt,wt2,
     &              20d0,430d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt lin',ptj1,wt,wt2,
     &              20d0,430d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 eta',etaj1,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt log',ptj2,wt,wt2,
     &              20d0,260d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt lin',ptj2,wt,wt2,
     &              20d0,260d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 eta',etaj2,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 3 pt log',ptj3,wt,wt2,
     &              20d0,170d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'Jet 3 pt lin',ptj3,wt,wt2,
     &              20d0,170d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 3 eta',etaj3,wt,wt2,
     &              -5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'(jet 1,jet 2) invariant mass',mjj,wt,wt2,
     &              0d0,425d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'delr',delr,wt,wt2,
     &              0.35d0,4.85d0,0.15d0,'lin')
      n=n+1
      call bookplot(n,tag,'deleta',deleta,wt,wt2,
     &              0d0,4d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'etastar',etastar,wt,wt2,
     &              0d0,3d0,0.6d0,'lin')
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
      
