      subroutine nplotter_Vgamma(p,wt,wt2,switch,nd)
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
      include 'frag.f'
      include 'phot_dip.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      double precision p(mxpart,4),wt,wt2
      double precision etarap,pt,r
      double precision y3,y4,y5,y6,pt3,pt4,pt5,pt6,re5,rea5,re6,rea6
      double precision yWgam,pt34,pttwo
      double precision yraptwo,m34,m345,yellgam,dot
      integer switch,n,nplotmax,nd 
      character*4 tag
      integer j,ilep,igam,inu
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
        pt3=1d7
        pt4=1d7
        pt34=1d7
c---Initialise photon 
        y5=1d3
        pt5=1d7
        yWgam=1d3
        yellgam=1d8
        m34=-1d0
        m345=-1d0
c----Initialise jet values will not pass cuts in there is an NLO jet
        y6=1d3
        pt6=1d7
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        re5=1d3
        rea5=1d3
        re6=1d3
        rea6=1d3
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
c--- If nproc=290, plot e^+(4). If nproc=295, plot e^-(3).
c--- For nproc=300 plot e^+(4) and e^-(3)
      if(nproc .eq. 290) then
         y4=etarap(4,p)
         pt4=pt(4,p)
      elseif(nproc .eq. 295) then 
         y3=etarap(3,p)
         pt3=pt(3,p)
      elseif(nproc .ge. 300) then 
         if(pt(3,p).gt.pt(4,p)) then 
            y3=etarap(3,p)
            pt3=pt(3,p)
            y4=etarap(4,p)
            pt4=pt(4,p)
         else
            y3=etarap(4,p)
            pt3=pt(4,p)
            y4=etarap(3,p)
            pt4=pt(3,p)
         endif
      endif
      
      pt34=pttwo(3,4,p)
      m34=dsqrt(max(2d0*dot(p,3,4),0d0))
      pt5=pt(5,p)
      y5=etarap(5,p)
      re5=-1d0
      rea5=-1d0
      if(nproc .eq. 290) then
         re5=R(p,4,5)
      elseif(nproc .eq. 295) then 
         re5=R(p,3,5)
      elseif((nproc .eq. 300) .or. (nproc .eq. 302)) then
         re5=R(p,3,5)
         rea5=R(p,4,5)
      endif
      
!---- Radiation Zero plot for all processes 
      yWgam=y5-yraptwo(3,4,p) 
      yellgam=1d8
      if(nproc.eq.290) then 
         yellgam=y5-y4     
      elseif(nproc.eq.295) then 
         yellgam=y5-y3
      endif

      if(jets .gt. 0) then
         pt6=pt(6,p)
         y6=etarap(6,p)
         re6=-1d0
         rea6=-1d0
         if(nproc .eq. 290) then
            re6=R(p,4,6)
         elseif(nproc .eq. 295) then 
            re6=R(p,3,6)
         elseif((nproc .eq. 300) .or. (nproc .eq. 302)) then
            re6=R(p,3,6)
            rea6=R(p,4,6)
         endif
      else ! put out of range of plotting
         pt6=1d7
         y6=1d7
         re6=1d7
         rea6=1d7
      endif


      if(nproc.lt.300) then 
c---  transverse mass of (e-gam,nu) system for Wgamma
        if ((nproc .eq. 290) .or. (nproc .eq. 292)) then
          inu=3
          ilep=4
        else
          inu=4
          ilep=3
        endif
        igam=5
        m345=(p(ilep,4)+p(igam,4))**2-(p(ilep,1)+p(igam,1))**2
     &      -(p(ilep,2)+p(igam,2))**2-(p(ilep,3)+p(igam,3))**2
        m345=m345+(p(ilep,1)+p(igam,1))**2
     &           +(p(ilep,2)+p(igam,2))**2
        m345=dsqrt(max(m345,0d0))+dsqrt(p(inu,1)**2+p(inu,2)**2)
        m345=m345**2
        do j=1,2 
           m345=m345-(p(3,j)+p(4,j)+p(5,j))**2
        enddo
        m345=dsqrt(max(m345,0d0)) 

c---  invariant mass of (Z,gam) system for Zgamma
      else 
         m345=(p(3,4)+p(4,4)+p(5,4))**2
         do j=1,3
           m345=m345-(p(3,j)+p(4,j)+p(5,j))**2 
         enddo        
         m345=dsqrt(max(m345,0d0))
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

      if((nproc .eq. 290).or.(nproc.ge.300)) then
         call bookplot(n,tag,'y4',y4,wt,wt2,-5d0,5d0,0.2d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt4',pt4,wt,wt2,0d0,100d0,2d0,'lin')
         n=n+1
      endif
      if((nproc.eq.295).or.(nproc.ge.300)) then
         call bookplot(n,tag,'y3',y3,wt,wt2,-5d0,5d0,0.2d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt3',pt3,wt,wt2,0d0,100d0,2d0,'lin')
         n=n+1
      endif
      call bookplot(n,tag,'boson invariant mass',
     & m34,wt,wt2,0d0,200d0,5d0,'lin')
      n=n+1
      
      call bookplot(n,tag,'DeltaRe5',re5,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
      if(nproc.eq.300) then 
      call bookplot(n,tag,'DeltaRea5',rea5,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
      endif
      call bookplot(n,tag,'pt34',pt34,wt,wt2,0d0,500d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0d0,500d0,10d0,'log')
      n=n+1
      if(nproc.lt.300) then 
      call bookplot(n,tag,'transverse cluster mass, m(e-gam,nu)',
     & m345,wt,wt2,0d0,200d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'ydiff(ellgam)',yellgam,wt,wt2,-5d0,5d0,
     & 0.2d0,'lin')
      n=n+1
      else
      call bookplot(n,tag,'(Z,gam) invariant mass',
     & m345,wt,wt2,0d0,200d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'(Z,gam) invariant mass',
     & m345,wt,wt2,0d0,4000d0,50d0,'log')
      n=n+1
      endif
      call bookplot(n,tag,'ydiff(Vgam)',yWgam,wt,wt2,-5d0,5d0,
     & 0.2d0,'lin')
      n=n+1      
      call bookplot(n,tag,'pt6',pt6,wt,wt2,0d0,200d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'y6',y6,wt,wt2,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaRe6',re6,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
c      if(nproc .lt. 305) then 
      call bookplot(n,tag,'DeltaRea6',rea6,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
c      endif

  
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
