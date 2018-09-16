      subroutine nplotter_gamgam(p,wt,wt2,switch,nd)
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
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,r,yraptwo,pttwo
      real(dp):: r2,phi34
      real(dp):: y3,y4,y34,y5,pt3,pt4,pt5,r45,r35,r34,s34,m34,pt34
      real(dp) :: mgg,ptgaga,pth,pts,ygh,ygs,ygg,yave
      common/gaga_observables/mgg,ptgaga,pth,pts,ygh,ygs,ygg
!$omp threadprivate(/gaga_observables/)
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
      real(dp) :: nj,swap
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/)
      real(dp):: gaga_highm_ptgaga(0:3),gaga_highm_ptgh(0:3)
      real(dp):: gaga_highm_ptgs(0:3),gaga_highm_ygg(0:3)
      real(dp):: gaga_highm_yh(0:3),gaga_highm_ys(0:3)
      common/gaga_highmgg_obs/gaga_highm_ptgaga,gaga_highm_ptgh,
     &     gaga_highm_ptgs,gaga_highm_ygg,gaga_highm_yh,gaga_highm_ys
!$omp threadprivate(/gaga_highmgg_obs/)
      real(dp) tiny
      parameter(tiny=1e-7_dp)
      
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        y3=-1d3
        y4=-1d3
        pt3=-1d3
        pt34=-1d3
        pt4=-1d3
        phi34=-1d3
        y34=-1d3
c---Intiailise jet 
        y5=1d3
        pt5=1d3
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        r34=1d3
        r35=1d3
        r45=1d3
        jets=nqcdjets
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
c--- Photons order based on pt
      pt3=pt(3,p)
      pt4=pt(4,p)
      if(pt3 > pt4) then 
         y3=yrap(3,p)
         y4=yrap(4,p)
       else
          swap=pt3
          pt3=pt4
          pt4=swap
          y3=yrap(4,p)
          y4=yrap(3,p)
       endif

!y34
       yave=0.5_dp*(y3+y4)
       y34=yraptwo(3,4,p)
       pt34=pttwo(3,4,p)

      
       r2= (p(3,1)*p(4,1)+p(3,2)*p(4,2))
     &     /(pt3*pt4)
       if (r2 > +0.9999999_dp) r2=+1._dp
       if (r2 < -0.9999999_dp) r2=-1._dp
       phi34=acos(r2)/(pi)
       if(abs(phi34-1._dp).lt.tiny)  phi34=phi34-tiny
!       write(6,*) phi34
!      pause
      
! R_ij
      r34=R(p,3,4)

! m_gamma,gamma
      s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
  
      m34=sqrt(s34)

!      if(m34.gt.700_dp) then
!         nj=jets
!      else
!         nj=-1
!      endif
      nj=real(jets,dp)
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

      call bookplot(n,tag,'mgaga',m34,wt,wt2,150._dp,1630._dp,40._dp
     &,'log')
      n=n+1

      call bookplot(n,tag,'ptga (hard)',pt3,wt
     &     ,wt2,0._dp,400._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptga (soft)',pt4,wt
     &     ,wt2,0._dp,400._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'y (hard)',y4,wt
     &     ,wt2,-3._dp,3._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y (soft)',y4,wt
     &     ,wt2,-3._dp,3._dp,0.4_dp,'lin')
      n=n+1

      call bookplot(n,tag,'ptgaga ',pt34,wt
     &     ,wt2,0._dp,1000._dp,25._dp,'log')
      n=n+1
      call bookplot(n,tag,'ygaga',y34,wt
     &     ,wt2,-3._dp,3._dp,0.4_dp,'lin')
      n=n+1

!---- observables for CMS comparision
!=====mgaga
      
      call bookplot(n,tag,'cms mgg 1',m34,wt
     &     ,wt2,0._dp,40._dp,40._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 2',m34,wt
     &     ,wt2,40._dp,60._dp,20._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 3',m34,wt
     &     ,wt2,60._dp,70._dp,10._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 4',m34,wt
     &     ,wt2,70._dp,100._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 5',m34,wt
     &     ,wt2,100._dp,120._dp,10._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 6',m34,wt
     &     ,wt2,120._dp,150._dp,30._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 7',m34,wt
     &     ,wt2,150._dp,250._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 8',m34,wt
     &     ,wt2,250._dp,400._dp,150._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms mgg 9',m34,wt
     &   ,wt2,400._dp,800._dp,400._dp,'lin')
      n=n+1


!=====ptgaga
      
      call bookplot(n,tag,'cms ptgg 1',pt34,wt
     &     ,wt2,0._dp,6._dp,6._dp,'lin')
      n=n+1

      call bookplot(n,tag,'cms ptgg 2',pt34,wt
     &     ,wt2,6._dp,10._dp,4._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 3',pt34,wt
     &     ,wt2,10._dp,24._dp,2._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 4',pt34,wt
     &     ,wt2,24._dp,28._dp,4._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 5',pt34,wt
     &     ,wt2,28._dp,40._dp,6._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 6',pt34,wt
     &     ,wt2,40._dp,100._dp,10._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 7',pt34,wt
     &     ,wt2,100._dp,120._dp,20._dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms ptgg 8',pt34,wt
     &     ,wt2,120._dp,200._dp,20._dp,'lin')
      n=n+1

!===== phi gaga
      call bookplot(n,tag,'cms phi 1',phi34,wt
     &     ,wt2,0._dp,0.6_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms phi 2',phi34,wt
     &     ,wt2,0.6_dp,0.8_dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms phi 3',phi34,wt
     &     ,wt2,0.8_dp,0.88_dp,0.04_dp,'lin')
      n=n+1
      call bookplot(n,tag,'cms phi 4',phi34,wt
     &     ,wt2,0.88_dp,1.00_dp,0.02_dp,'lin')
      n=n+1
      
!-----observables for high mgg study

      call bookplot(n,tag,'ptgaga(inc)',gaga_highm_ptgaga(0),wt
     &   ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptgaga(120-130)',gaga_highm_ptgaga(1),wt
     &   ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
!      call bookplot(n,tag,'ptgaga(700-800)',gaga_highm_ptgaga(2),wt
!     &   ,wt2,0._dp,500._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'ptgaga(800+)',gaga_highm_ptgaga(3),wt
!     &   ,wt2,0._dp,500._dp,10._dp,'log')
!      n=n+1

      call bookplot(n,tag,'ptgah(inc)',gaga_highm_ptgh(0),wt
     &   ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptgah(120-130)',gaga_highm_ptgh(1),wt
     &   ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
!      call bookplot(n,tag,'ptgh(700-800)',gaga_highm_ptgh(2),wt
!     &   ,wt2,100._dp,600._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'ptgah(800+)',gaga_highm_ptgh(3),wt
!     &   ,wt2,200._dp,700._dp,10._dp,'log')
!      n=n+1


      call bookplot(n,tag,'ptgas(inc)',gaga_highm_ptgs(0),wt
     &     ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptgs(120-130)',gaga_highm_ptgs(1),wt
     &   ,wt2,0._dp,200._dp,4._dp,'log')
      n=n+1
!      call bookplot(n,tag,'ptgs(700-800)',gaga_highm_ptgs(2),wt
!     &   ,wt2,100._dp,600._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'ptgas(800+)',gaga_highm_ptgs(3),wt
!     &   ,wt2,200._dp,700._dp,10._dp,'log')
!      n=n+1

      call bookplot(n,tag,'ygaga(inc)',gaga_highm_ygg(0),wt
     &     ,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'ygaga(120-130)',gaga_highm_ygg(1),wt
     &     ,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1
!      call bookplot(n,tag,'ygaga(700-800)',gaga_highm_ygg(2),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygaga(800)',gaga_highm_ygg(3),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1


!      call bookplot(n,tag,'ygah(inc)',gaga_highm_yh(0),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygah(200-700)',gaga_highm_yh(1),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygah(700-800)',gaga_highm_yh(2),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygah(800)',gaga_highm_yh(3),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1!

!      call bookplot(n,tag,'ygas(inc)',gaga_highm_ys(0),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygas(200-700)',gaga_highm_ys(1),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygas(700-800)',gaga_highm_ys(2),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'ygas(800)',gaga_highm_ys(3),wt
!     &     ,wt2,-4._dp,4._dp,0.2_dp,'lin')
!      n=n+1
    
      call bookplot(n,tag,'|y(3)|',abs(y3),wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'|y(4)|',abs(y4),wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'|y(3)| and |y(4)|',abs(y3),half*wt,quarter*wt2,0._dp,5._dp,0.1_dp,'lin')
      call bookplot(n,tag,'|y(3)| and |y(4)|',abs(y4),half*wt,quarter*wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1

      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,500._dp,10._dp,'lin')
      n=n+1

      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,500._dp,10._dp,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-4._dp,4._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Abs(y34)',abs(y34),wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1

      
      call bookplot(n,tag,'(y3+y4)/2',yave,wt,wt2,
     & -6._dp,6._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,2._dp,202._dp,4._dp,'log')
      n=n+1

      call bookplot(n,tag,'pt3_2',pt3,wt,wt2,0._dp,80._dp,2.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,2._dp,202._dp,4._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt4_2',pt4,wt,wt2,0._dp,100._dp,2.5_dp,'lin')
      n=n+1
                
      call bookplot(n,tag,'y34',y34,wt,wt2,-6._dp,6._dp,0.5_dp,'lin')
      n=n+1
     
      call bookplot(n,tag,'DeltaR34',r34,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
     
!    m34 
      call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,300._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,500._dp,1000._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,'m34',m34,wt,wt2,300._dp,1000._dp,20._dp,'log')
      n=n+1
      call bookplot(n,tag,'nj',nj,wt,wt2,-0.5_dp,3.5_dp,1._dp,'lin')
      n=n+1

! Extra plots
!      call bookplot(n,tag,'m34',m34,wt,wt2,200._dp,1600._dp,46.667_dp,'log')
!      n=n+1
!      call bookplot(n,tag,'pth',pt3,wt,wt2,0._dp,1000._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'pts',pt4,wt,wt2,0._dp,1000._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,1000._dp,10._dp,'log')
!      n=n+1
!      call bookplot(n,tag,'phi34',phi34*pi,wt,wt2,0._dp,3.1415936_dp,0.0981748_dp,'lin')
!      n=n+1
      
      call bookplot(n,tag,'y5',y5,wt,wt2,-3._dp,3._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,200._dp,2._dp,'log')
      n=n+1
   
!      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
!      n=n+1
!      call bookplot(n,tag,'DeltaR45',r45,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
!      n=n+1
  
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
      
      

