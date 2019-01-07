cccccccccccccc GENERIC NPLOTTER FOR ZGAMJET ccccccccccccccccc

      subroutine nplotter_zgamjet(p,wt,wt2,switch,nd)
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
      include 'outputflags.f'
      include 'nproc.f'
c-----
      character*4 mypart
      common/mypart/mypart
c-----
      double precision p(mxpart,4),wt,wt2
      double precision r,etmiss
      double precision s34,m34
c-----
      double precision ptgam
      double precision m345,ptmiss,etvec(4)
      double precision Raj1,Raj2,Rajmin,Ralm

c-----
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
        m34=0D0
        m345=0D0
        ptgam=0D0
        Raj1=0D0
        Raj2=0D0
        Rajmin=0D0
        ptmiss=0D0
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
c-1---m(l,l)
      s34=2d0*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
      m34=dsqrt(s34)
c-2---m(l,l,gamma)
      m345=dsqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     .          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
c-3----pT(photon)
       ptgam = dsqrt(p(5,1)**2+p(5,2)**2)
c-5-6--R(gam,jet)
       Raj1 = R(p,5,6)
       Raj2 = R(p,5,7)
       if (mypart.eq.'tota') then
          Rajmin = min(Raj1,Raj2)
       else
          Rajmin = Raj1
       endif
       if (nproc.eq.302) then
c------only electron case
c-7-------lepton angular difference
c          deltall = dabs(deltaPHI(p,3,4))
c-8-------R(jet,lepton)
c          Rlmj1 = R(p,3,6)
c          Rlmj2 = R(p,3,7)
c          if (mypart.eq.'tota') then
c             Rlmjmin = min(Rlmj1,Rlmj2)
c          else
c             Rlmjmin = Rlmj1
c          endif
c-9-------R(jet,antilepton)
c          Rlpj1 = R(p,4,6)
c          Rlpj2 = R(p,4,7)
c          if (mypart.eq.'tota') then
c             Rlpjmin = min(Rlpj1,Rlpj2)
c          else
c             Rlpjmin = Rlpj1
c          endif
c-10------R(gamma,lepton)
          Ralm = R(p,5,3)
c-11------R(gamma,antilepton)
c          Ralp = R(p,5,4)
c-12------delta eta(gamma,lepton)
c          detalma = dabs(etarap(3,p)-etarap(5,p))
c-13------delta eta(gamma,antilepton)
c          detalpa = dabs(etarap(4,p)-etarap(5,p))
       else if (nproc.eq.307) then
c------only neutrino case
c---------missing pT
          ptmiss=etmiss(p,etvec)
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

!-1---m34 
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
!-2---m345
      call bookplot(n,tag,'m(l,l,gam)',m345,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
!-3---photon pT
      call bookplot(n,tag,'pT(gam)',ptgam,wt,wt2,15d0,400d0,4d0,'lin')
      n=n+1
!-3b--photon pT
      call bookplot(n,tag,'pT(gam)',ptgam,wt,wt2,100d0,400d0,4d0,'lin')
      n=n+1
!-4---pT-jet max
c      call bookplot(n,tag,'pT(j)max',ptjmax,wt,wt2,0d0,400d0,10d0,'lin')
c      n=n+1
!-5a---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,10d0,400d0,10d0,'lin')
c      n=n+1
!-5b---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,15d0,400d0,10d0,'lin')
c      n=n+1
!-5c---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,20d0,400d0,10d0,'lin')
c      n=n+1
!-6---min R(gam,jet) [1]
      call bookplot(n,tag,'Rajmin',Rajmin,wt,wt2,0d0,4d0,0.08d0,'lin')
      n=n+1
!-7---min R(gam,jet) [2]
      call bookplot(n,tag,'Rajmin',Rajmin,wt,wt2,0d0,1.4d0,0.04d0,'lin')
      n=n+1
c-----
      if (nproc.eq.302) then
!-8---delta phi(l,l)
c      call bookplot(n,tag,'dphi_ll',deltall,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
!-9---min R(j,l-)
c      call bookplot(n,tag,'min_Rjlm',Rlmjmin,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
!-10--min R(j,l+)
c      call bookplot(n,tag,'min_Rjlp',Rlpjmin,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
!-11--R(gam,l-)
      call bookplot(n,tag,'min_Ralm',Ralm,wt,wt2,0d0,5d0,0.1d0,'lin')
      n=n+1
!-11b--R(gam,l-)
      call bookplot(n,tag,'min_Ralm',Ralm,wt,wt2,1.5d0,5d0,0.1d0,'lin')
      n=n+1
!-12--R(gam,l+)
c      call bookplot(n,tag,'min_Ralp',Ralp,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
!-13--|eta(gam)-eta(l-)|
c      call bookplot(n,tag,'deta_alm',detalma,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
!-14--|eta(gam)-eta(l+)|
c      call bookplot(n,tag,'deta_alp',detalpa,wt,wt2,0d0,5d0,0.1d0,'lin')
c      n=n+1
c-----
      elseif (nproc.eq.307) then      
!-7---missing transverse momentum
      call bookplot(n,tag,'pT(miss)',ptmiss,wt,wt2,25d0,400d0,4d0,'lin')
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
      
      
************************************************************************
*      delta PHI                                                       *
************************************************************************
      double precision function deltaPHI(p,i,j)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),r2
      integer i,j
c-----
      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))
     .     /dsqrt((p(i,1)**2+p(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      deltaPHI=dacos(r2)
c-----
      return
      end


cccccccccccccc NPLOTTER FOR ZGAMJET: Xsec table w PT(Gam) bin  ccccccccccccccccc
      subroutine nplotter_zgamjet_ptgam(p,wt,wt2,switch,nd)
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'outputflags.f'
c-----
      double precision p(mxpart,4),wt,wt2,ptgam
c-----
      integer switch,n,nplotmax
      character*4 tag
      integer nd
      logical first
      common/nplotmax/nplotmax
      data first/.true./
      save first
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        ptgam=0D0
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
c------pT(photon)
       ptgam = dsqrt(p(5,1)**2+p(5,2)**2)

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
      n=1              

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

!-----photon pT
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,15d0,10000d0,100d0,'lin')
      n=n+1
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,30d0,10000d0,100d0,'lin')
      n=n+1
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,60d0,10000d0,100d0,'lin')
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
 

cccccccccccccc NPLOTTER FOR ZGAMJET: Exclusive Dist.  ccccccccccccccccc
      subroutine nplotter_zgam_1j(p,wt,wt2,switch,nd)
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'outputflags.f'
      double precision p(mxpart,4),wt,wt2
      double precision s34,m34,ptgam,m345
c-----
      integer switch,n,nplotmax
      character*4 tag
      integer nd
      logical first
      common/nplotmax/nplotmax
      data first/.true./
      save first
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        m34=0D0
        m345=0D0
        ptgam=0D0
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
c------pT(photon)
       ptgam = dsqrt(p(5,1)**2+p(5,2)**2)
c-----m(l,l)
      s34=2d0*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
      m34=dsqrt(s34)
c-----m(l,l,gamma)
      m345=dsqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)

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
      n=1              

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

!-----photon pT
      call bookplot(n,tag,'pT(ga15)',ptgam,wt,wt2,15d0,300d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'pT(ga60)',ptgam,wt,wt2,60d0,300d0,10d0,'lin')
      n=n+1
!-----mll
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0d0,250d0,5d0,'lin')
      n=n+1
!-----mllgam
      call bookplot(n,tag,'m(l,l,gam)',m345,wt,wt2,0d0,250d0,5d0,'lin')
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
 
