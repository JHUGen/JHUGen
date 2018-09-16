cccccccccccccc GENERIC NPLOTTER FOR ZGAMJET ccccccccccccccccc

      subroutine nplotter_zgamjet(p,wt,wt2,switch,nd)
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
!-- nd determines whether we are analysing a photon dipole and hence have
!---to rescale accordingly
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      include 'nproc.f'
      include 'kpart.f'
c-----
      integer mykpart
      common/mykpart/mykpart
c-----
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: r,etmiss
      real(dp):: s34,m34
c-----
      real(dp):: ptgam
      real(dp):: m345,ptmiss,etvec(4)
      real(dp):: Raj1,Raj2,Rajmin,Ralm

c-----
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
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
        m34=0._dp
        m345=0._dp
        ptgam=0._dp
        Raj1=0._dp
        Raj2=0._dp
        Rajmin=0._dp
        ptmiss=0._dp
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
c-1---m(l,l)
      s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
      m34=sqrt(s34)
c-2---m(l,l,gamma)
      m345=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
c-3----pT(photon)
       ptgam = sqrt(p(5,1)**2+p(5,2)**2)
c-5-6--R(gam,jet)
       Raj1 = R(p,5,6)
       Raj2 = R(p,5,7)
       if (mykpart==ktota) then
          Rajmin = min(Raj1,Raj2)
       else
          Rajmin = Raj1
       endif
       if (nproc==302) then
c------only electron case
c-7-------lepton angular difference
c          deltall = abs(deltaPHI(p,3,4))
c-8-------R(jet,lepton)
c          Rlmj1 = R(p,3,6)
c          Rlmj2 = R(p,3,7)
c          if (mykpart==ktota) then
c             Rlmjmin = min(Rlmj1,Rlmj2)
c          else
c             Rlmjmin = Rlmj1
c          endif
c-9-------R(jet,antilepton)
c          Rlpj1 = R(p,4,6)
c          Rlpj2 = R(p,4,7)
c          if (mykpart==ktota) then
c             Rlpjmin = min(Rlpj1,Rlpj2)
c          else
c             Rlpjmin = Rlpj1
c          endif
c-10------R(gamma,lepton)
          Ralm = R(p,5,3)
c-11------R(gamma,antilepton)
c          Ralp = R(p,5,4)
c-12------delta eta(gamma,lepton)
c          detalma = abs(etarap(3,p)-etarap(5,p))
c-13------delta eta(gamma,antilepton)
c          detalpa = abs(etarap(4,p)-etarap(5,p))
       else if (nproc==307) then
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

!-1---m34 
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0._dp,250._dp,5._dp,'lin')
      n=n+1
!-2---m345
      call bookplot(n,tag,'m(l,l,gam)',m345,wt,wt2,0._dp,250._dp,5._dp,'lin')
      n=n+1
!-3---photon pT
      call bookplot(n,tag,'pT(gam)',ptgam,wt,wt2,15._dp,400._dp,4._dp,'lin')
      n=n+1
!-3b--photon pT
      call bookplot(n,tag,'pT(gam)',ptgam,wt,wt2,100._dp,400._dp,4._dp,'lin')
      n=n+1
!-4---pT-jet max
c      call bookplot(n,tag,'pT(j)max',ptjmax,wt,wt2,0._dp,400._dp,10._dp,'lin')
c      n=n+1
!-5a---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,10._dp,400._dp,10._dp,'lin')
c      n=n+1
!-5b---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,15._dp,400._dp,10._dp,'lin')
c      n=n+1
!-5c---pT-jet min
c      call bookplot(n,tag,'pTj min',ptjmin,wt,wt2,20._dp,400._dp,10._dp,'lin')
c      n=n+1
!-6---min R(gam,jet) [1]
      call bookplot(n,tag,'Rajmin',Rajmin,wt,wt2,0._dp,4._dp,0.08_dp,'lin')
      n=n+1
!-7---min R(gam,jet) [2]
      call bookplot(n,tag,'Rajmin',Rajmin,wt,wt2,0._dp,1.4_dp,0.04_dp,'lin')
      n=n+1
c-----
      if (nproc==302) then
!-8---delta phi(l,l)
c      call bookplot(n,tag,'dphi_ll',deltall,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
!-9---min R(j,l-)
c      call bookplot(n,tag,'min_Rjlm',Rlmjmin,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
!-10--min R(j,l+)
c      call bookplot(n,tag,'min_Rjlp',Rlpjmin,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
!-11--R(gam,l-)
      call bookplot(n,tag,'min_Ralm',Ralm,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
!-11b--R(gam,l-)
      call bookplot(n,tag,'min_Ralm',Ralm,wt,wt2,1.5_dp,5._dp,0.1_dp,'lin')
      n=n+1
!-12--R(gam,l+)
c      call bookplot(n,tag,'min_Ralp',Ralp,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
!-13--|eta(gam)-eta(l-)|
c      call bookplot(n,tag,'deta_alm',detalma,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
!-14--|eta(gam)-eta(l+)|
c      call bookplot(n,tag,'deta_alp',detalpa,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
c-----
      elseif (nproc==307) then      
!-7---missing transverse momentum
      call bookplot(n,tag,'pT(miss)',ptmiss,wt,wt2,25._dp,400._dp,4._dp,'lin')
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
      function deltaPHI(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: deltaPHI
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),r2
      integer:: i,j
c-----
      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))
     &     /sqrt((p(i,1)**2+p(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      deltaPHI=acos(r2)
c-----
      return
      end


cccccccccccccc NPLOTTER FOR ZGAMJET: Xsec table w PT(Gam) bin  ccccccccccccccccc
      subroutine nplotter_zgamjet_ptgam(p,wt,wt2,switch,nd)
      implicit none
      include 'types.f'
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
c-----
      real(dp):: p(mxpart,4),wt,wt2,ptgam
c-----
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
      logical:: first
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
        tag=tagbook
        ptgam=0._dp
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
c------pT(photon)
       ptgam = sqrt(p(5,1)**2+p(5,2)**2)

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
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,15._dp,10000._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,30._dp,10000._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pT(A)',ptgam,wt,wt2,60._dp,10000._dp,100._dp,'lin')
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
      include 'types.f'
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: s34,m34,ptgam,m345
c-----
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
      logical:: first
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
        tag=tagbook
        m34=0._dp
        m345=0._dp
        ptgam=0._dp
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
c------pT(photon)
       ptgam = sqrt(p(5,1)**2+p(5,2)**2)
c-----m(l,l)
      s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
      m34=sqrt(s34)
c-----m(l,l,gamma)
      m345=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
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
        call bookfill(tag,p,wt/real(itmx,dp))  
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
      call bookplot(n,tag,'pT(ga15)',ptgam,wt,wt2,15._dp,300._dp,10._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pT(ga60)',ptgam,wt,wt2,60._dp,300._dp,10._dp,'lin')
      n=n+1
!-----mll
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0._dp,250._dp,5._dp,'lin')
      n=n+1
!-----mllgam
      call bookplot(n,tag,'m(l,l,gam)',m345,wt,wt2,0._dp,250._dp,5._dp,'lin')
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
 
