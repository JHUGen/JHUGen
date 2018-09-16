      subroutine nplotter_Zbbbar(p,wt,wt2,switch,nd)
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
c---   nd:    dipole number if appropriate      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'jetlabel.f'
      include 'outputflags.f'
      include 'nproc.f'
      integer:: inu,ibq,iba,ilt,nd
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,r,yraptwo,pttwo
      real(dp):: ylept1,ylept2,ptlept1,ptlept2,yw,ptw,mw,ybq,ptbq,yba,ptba,
     & ylt,ptlt,yj1,ptj1,yj2,ptj2,ybmax,yj3,ptj3,mbb,mjj
     & ,dphibqnu,dphibanu,r2,binedges(30),tag2b,ptbmin,maxyb
      real(dp):: etaw, etabq, etaba, etaraptwo
      integer:: switch,n,nplotmax
      integer tag
      common/nplotmax/nplotmax
      logical, save::first=.true.
ccccc!$omp threadprivate(first,/nplotmax/)

  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ylept1=99._dp
      ylept2=99._dp
      yw=99._dp
      ybq=99._dp
      yba=99._dp
      ylt=99._dp
      yj1=99._dp
      yj2=99._dp
      yj3=99._dp
      ptlept1=-1._dp
      ptlept2=-1._dp
      ptw=-1._dp
      ptbq=-1._dp
      ptba=-1._dp
      ptlt=-1._dp
      ptj1=-1._dp
      ptj2=-1._dp
      ptj3=-1._dp
      mw=-1._dp
      mbb=-1._dp
      mjj=-1._dp
      dphibqnu=99._dp
      dphibanu=99._dp
      etaw=99._dp
      etabq=99._dp
      etaba=99._dp
c      ybmax=99._dp

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
*     Relevant processes are:                                          *
*                                                                      *
* 20  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6) [massive]   *
* 25  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6) [massive]  *
*                                                                      *
************************************************************************
      
      ylept1=yrap(3,p)
      ptlept1=pt(3,p)
      ylept2=yrap(4,p)
      ptlept2=pt(4,p)
      
      yw=yraptwo(3,4,p)
      etaw=etaraptwo(3,4,p)
      ptw=pttwo(3,4,p)
      mw=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)


c--- find the locations of the b,anti-b and light (or gluon) quarks      
      ibq=-1
      iba=-1
      ilt=-1
      if     (jetlabel(1) == 'bq') then
        ibq=5
      elseif ((jetlabel(2) == 'bq') .and. (jets >= 2)) then
        ibq=6
      elseif ((jetlabel(3) == 'bq') .and. (jets >= 3)) then
        ibq=7
      endif
      
      if     (jetlabel(1) == 'ba') then
        iba=5
      elseif ((jetlabel(2) == 'ba') .and. (jets >= 2)) then
        iba=6
      elseif ((jetlabel(3) == 'ba') .and. (jets >= 3)) then
        iba=7
      endif
      
      if     (jetlabel(1) == 'pp') then
        ilt=5
      elseif ((jetlabel(2) == 'pp') .and. (jets >= 2)) then
        ilt=6
      elseif ((jetlabel(3) == 'pp') .and. (jets >= 3)) then
        ilt=7
      endif
      
c--- this case only occurs for processes 421, 426
      if     (jetlabel(1) == 'bb') then
        ibq=5
      elseif ((jetlabel(2) == 'bb') .and. (jets >= 2)) then
        ibq=6
      endif
      
      if (ibq > 0) then
        ybq=yrap(ibq,p)
        ptbq=pt(ibq,p)
      endif
      if (iba > 0) then
        yba=yrap(iba,p)
        ptba=pt(iba,p)
      endif
      if (ilt > 0) then
        ylt=yrap(ilt,p)
        ptlt=pt(ilt,p)
      endif
       
      yj1=yrap(5,p)
      ptj1=pt(5,p)
      if (jets >= 2) then
        yj2=yrap(6,p)
        ptj2=pt(6,p)
      endif
      if (jets >= 3) then
        yj3=yrap(6,p)
        ptj3=pt(6,p)
      endif

      if (jets > 1) then
        mjj=sqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &           -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
      endif
      if ((ibq > 0) .and. (iba > 0)) then
        mbb=sqrt((p(ibq,4)+p(iba,4))**2-(p(ibq,1)+p(iba,1))**2
     &           -(p(ibq,2)+p(iba,2))**2-(p(ibq,3)+p(iba,3))**2)
      endif

c---  Add deltaphi histogram for bq and ba
      if (ibq > 0) then
         r2 = (p(ibq,1)*p(inu,1)+p(ibq,2)*p(inu,2))/
     & (sqrt((p(ibq,1)**2+p(ibq,2)**2)*(p(inu,1)**2+p(inu,2)**2)))
         if (r2 > +0.9999999_dp) r2=+1._dp
         if (r2 < -0.9999999_dp) r2=-1._dp
         dphibqnu = acos(r2)
      endif

      if (iba > 0) then
         r2 = (p(iba,1)*p(inu,1)+p(iba,2)*p(inu,2))/
     & (sqrt((p(iba,1)**2+p(iba,2)**2)*(p(inu,1)**2+p(inu,2)**2)))
         if (r2 > +0.9999999_dp) r2=+1._dp
         if (r2 < -0.9999999_dp) r2=-1._dp
         dphibanu = acos(r2)
      endif


c--   pt and rapidity distributions for the highest pt b-quark jet
      if ((iba < 0) .and. (ibq > 0)) then
c         ptbmax = pt(ibq, p)
         ybmax = yrap(ibq, p)
      elseif ((iba > 0) .and. (ibq < 0)) then
c         ptbmax = pt(iba, p)
         ybmax = yrap(iba, p)
      elseif ((iba > 0) .and. (ibq > 0) .and. 
     &    (pt(iba, p) > pt(ibq, p)) ) then
c         ptbmax = pt(iba, p)
         ybmax = yrap(iba, p)
      elseif ((iba > 0) .and. (ibq > 0) .and.
     &    (pt(iba,p) <= pt(ibq, p)) ) then
c         ptbmax = pt(ibq, p)
         ybmax = yrap(ibq, p)
      endif

         

      
c      write(6,*) 'In Wbbmas, switch=',switch
c      write(6,*) 'In Wbbmas, jets=',jets
c      write(6,*) 'In Wbbmas, jetlabel=',jetlabel
c      write(6,*)
      
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

       ptbmin=min(ptbq,ptba)

      call bookplot(n,tag,'softest b-jet pt',ptbmin,wt,wt2,
     &              0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'+INTEGRAL+ softest b-jet pt',ptbmin,wt,wt2,
     &              0._dp,500._dp,10._dp,'log')
      n=n+1

      call bookplot(n,tag,'Lept 1 rapidity',ylept1,wt,wt2,
     &              -10._dp,10._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Lept 1 pt',ptlept1,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      n=n+1

      call bookplot(n,tag,'Lept 2 rapidity',ylept2,wt,wt2,
     &              -10._dp,10._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Lept 2 pt',ptlept2,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      n=n+1

      call bookplot(n,tag,'Z boson rapidity',yw,wt,wt2,
     &              -10._dp,10._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Z boson pt',ptw,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'Z boson invariant mass',mw,wt,wt2,
     &              0._dp,200._dp,5._dp,'lin')
      n=n+1

      maxyb=max(abs(ybq),abs(yba))
      call bookplot(n,tag,'max abs. b-quark jet rapidity',maxyb,wt,wt2,
     &              0._dp,10._dp,0.2_dp,'lin')
      n=n+1

      call bookplot(n,tag,'b-quark jet rapidity',ybq,wt,wt2,
     &              -10._dp,10._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'b-quark jet pt',ptbq,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'anti-b-quark jet rapidity',yba,wt,wt2,
     &              -10._dp,10._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'anti-b-quark jet pt',ptba,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'(b,anti-b) jet invariant mass',mbb,wt,wt2,
     &              0._dp,2000._dp,20._dp,'log')
      n=n+1

      call bookplot(n,tag,'Jet 1 rapidity',yj1,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 2 rapidity',yj2,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 2 pt',ptj2,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'(jet 1,jet 2) invariant mass',mjj,wt,wt2,
     &              0._dp,800._dp,10._dp,'log')
      n=n+1

      call bookplot(n,tag,'light jet rapidity',ylt,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'light jet pt',ptlt,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 3 rapidity',yj3,wt,wt2,
     &              -5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Jet 3 pt',ptj3,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1

c---  Add histograms for dphi_bqnu and dphi_banu

      call bookplot(n,tag,'deltaphi(nu,b-quark-jet)', dphibqnu, 
     & wt, wt2, 0._dp, pi+0.1_dp, 0.1_dp, 'lin')
      n=n+1

      call bookplot(n,tag,'deltaphi(nu,anti-b-quark-jet)', dphibanu, 
     & wt, wt2, 0._dp, pi+0.1_dp, 0.1_dp, 'lin')
      n=n+1

c---  Add pseudorapidity histograms

      call bookplot(n,tag,'Z boson pseudo-rapidity',etaw,wt,wt2,
     &              -5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'b-quark jet pseudo-rapidity',etabq,wt,wt2,
     &              -5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'anti-b-quark jet pseudo-rapidity',etaba,
     & wt,wt2, -5._dp,5._dp,0.1_dp,'lin')
      n=n+1

c---   pt and rapidity for the highest pt b-quark jet

      call bookplot(n,tag,'b-quark harder jet rapidity',ybq,wt,wt2,
     &              -5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'b-quark harder jet pt',ptbq,wt,wt2,
     &              0._dp,200._dp,5._dp,'log')
      n=n+1

      call bookplot(n,tag,'b-jet pt (b and b-bar binned)',ptbq,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
      call bookplot(n,tag,'b-jet pt (b and b-bar binned)',ptba,wt,wt2,
     &              0._dp,500._dp,5._dp,'log')
c-- Irregular bins: note that binedges must have 30 entries
c      if (first) then
c       binedges= (/ 20._dp, 30._dp, 40._dp, 50._dp, 70._dp, 150._dp, 500._dp,
c     &             (0._dp,j=1,23) /)
c       call initirregbins(n,binedges)
c      endif
      n=n+1

      tag2b=0._dp
      if ((ibq > 0) .and. (iba > 0)) tag2b=1._dp
      call bookplot(n,tag,'b-jet pt (2b events only)',ptbq,
     &              wt*tag2b,wt2*tag2b**2,0._dp,500._dp,5._dp,'log')
      call bookplot(n,tag,'b-jet pt (2b events only)',ptba,
     &              wt*tag2b,wt2*tag2b**2,0._dp,500._dp,5._dp,'log')
c-- Irregular bins: note that binedges must have 30 entries
c      if (first) then
c       binedges= (/ 20._dp, 30._dp, 40._dp, 50._dp, 70._dp, 150._dp, 500._dp,
c     &             (0._dp,j=1,23) /)
c       call initirregbins(n,binedges)
c      endif
      n=n+1
      
c      call bookplot(n,tag,'DeltaRe5',re5,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'y5',y5,wt,wt2,-yjet,yjet,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,ptjet,2._dp,'lin')
c      n=n+1

  
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
      
