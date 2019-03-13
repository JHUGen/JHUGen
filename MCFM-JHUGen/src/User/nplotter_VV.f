      subroutine nplotter_VV(p,wt,wt2,switch,nd)
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
      include 'masses.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      double precision p(mxpart,4),wt,wt2
      double precision etarap,pt
      double precision y4,y5,pt3,pt4,pt5,pt6
      double precision pt34,pttwo
      double precision pt45,pt56,m45,mt45,mtrans,Etll,MET
      double precision etvec(4),etmiss,r2,delphi,m34,m56,m3456
      double precision mthl,mthu
      integer switch,n,nplotmax,nd 
      character*4 tag
      integer i
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/,y4,pt3,pt4,y5)

!=====info for "peak" transverse mass plot
      mthu=hmass*1d0
      mthl=0.75d0*hmass
      
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        y4=1d3
        pt3=1d7
        pt4=1d7
        pt34=1d7
c---Intiailise photon 
        y5=1d3
        pt5=1d7
        mt45=1d7
        m45=1d7
        pt45=1d7
        pt34=1d7
        pt56=1d7
c----Initialise jet values will not pass cuts in there is an NLO jet
        pt6=1d7
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

!     121 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [top, bottom loops, exact]' 'L'
!     122 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [above + interf. with gg->WW]' 'L'
      
         pt3=pt(3,p)
         pt4=pt(4,p) 
         pt5=pt(5,p)
         pt6=pt(6,p) 
         y4=etarap(4,p)
         y5=etarap(5,p) 
         pt45=pttwo(4,5,p)
         pt34=pttwo(3,4,p) 
         pt56=pttwo(5,6,p)
         r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .        /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
         if (r2 .gt. +0.9999999D0) r2=+1D0
         if (r2 .lt. -0.9999999D0) r2=-1D0
         delphi=dacos(r2)
!--- mll cut 
         m45=(p(4,4)+p(5,4))**2 
         do i=1,3
           m45=m45-(p(4,i)+p(5,i))**2
         enddo
         m45=dsqrt(m45)
         m34=(p(3,4)+p(4,4))**2 
         m3456=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
         do i=1,3
           m34=m34-(p(3,i)+p(4,i))**2
           m3456=m3456-(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2
        enddo
        m34=dsqrt(m34)
        m3456=dsqrt(m3456)
        m56=(p(5,4)+p(6,4))**2 
        do i=1,3
           m56=m56-(p(5,i)+p(6,i))**2
        enddo
        m56=dsqrt(m56)

      mt45=0d0 
      mt45=(dsqrt(dsqrt(pttwo(4,5,p)**2+m45**2)+etmiss(p,etvec))**2)

      etvec(:)=p(3,:)+p(6,:)
c--- transverse mass
      Etll=dsqrt(max(0d0,(p(4,4)+p(5,4))**2-(p(4,3)+p(5,3))**2))
      MET=dsqrt(max(0d0,etvec(1)**2+etvec(2)**2))
      mtrans=MET+Etll
     
c      Etll=dsqrt(
c     . +(p(4,4)+p(5,4))**2
c     . -(p(4,3)+p(5,3))**2)
c      MET=dsqrt(
c     . +etvec(1)**2+etvec(2)**2
c     . +(p(4,4)+p(5,4))**2
c     . -(p(4,1)+p(5,1))**2
c     . -(p(4,2)+p(5,2))**2
c     . -(p(4,3)+p(5,3))**2)
c      mtrans=dsqrt((MET+Etll)**2
c     . -(p(4,1)+p(5,1)+etvec(1))**2
c     . -(p(4,2)+p(5,2)+etvec(2))**2)

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
 
      call bookplot(n,tag,'pt_nu_1',pt3,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_e_1',pt4,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y_e_1',y4,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_W_1',pt34,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_e_2',pt5,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'y_e_2',y5,wt,wt2,-4d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_nu_2',pt6,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_W_2',pt56,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt_ll',pt45,wt,wt2,0d0,100d0,2.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'mll',m45,wt,wt2,0d0,100d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'mll',m34,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'mll_2',m56,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
c      call bookplot(n,tag,'m4l',m3456,wt,wt2,hmass-0.5d0,hmass+0.5d0,
c     & 0.05d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m4l',m3456,wt,wt2,0d0,2000d0,20d0
c     &,'log')
c      n=n+1
      call bookplot(n,tag,'mt',mt45,wt,wt2,0d0,1000d0,20d0,'lin')
      n=n+1
      call bookplot(n,tag,'delphi',delphi,wt,wt2,0d0,3.14d0,0.1d0,'lin')
      n=n+1
  
      call bookplot(n,tag,'transverse mass',
     & mtrans,wt,wt2,10d0,510d0,10d0,'log')
      n=n+1

c--- Plots of mtrans in specific regions
      call bookplot(n,tag,'10 < m(trans) < 2010',
     & mtrans,wt,wt2,10d0,2010d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < m(trans) < 2010',
     & mtrans,wt,wt2,130d0,2010d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < m(trans) < 2020',
     & mtrans,wt,wt2,300d0,2020d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(trans) < 130',
     & mtrans,wt,wt2,10d0,130d0,5d0,'lin')
      n=n+1

c--- Plots of mtrans in specific regions
      call bookplot(n,tag,'.75 m_H < m(trans) < m_H',
     & mtrans,wt,wt2,mthl,mthu,1d0,'log')
      n=n+1
      
      
c--- Plots of m(3456) in specific regions
      call bookplot(n,tag,'10 < m(3456) < 2010',
     & m3456,wt,wt2,10d0,2010d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < m(3456) < 2010',
     & m3456,wt,wt2,130d0,2010d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < m(3456) < 2020',
     & m3456,wt,wt2,300d0,2020d0,20d0,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(3456) < 130',
     & m3456,wt,wt2,10d0,130d0,5d0,'lin')
      n=n+1
      
      call bookplot(n,tag,'pt(W)',
     & pt34,wt,wt2,0d0,2d0,0.02d0,'lin')
      n=n+1
      
      call bookplot(n,tag,'+INTEGRAL+ pt(W)',
     & pt34,wt,wt2,0d0,10d0,0.1d0,'lin')
      n=n+1
      
      call bookplot(n,tag,'50 < m(3456) < 250',
     & m3456,wt,wt2,50d0,250d0,2d0,'log')
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
      
c
      return 
      end
