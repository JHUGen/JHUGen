      subroutine nplotter_qqZZqq(p,wt,wt2,switch)
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
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      include 'masses.f'
      include 'removebr.f'
      real(dp):: p(mxpart,4),wt,wt2,m3456,pt34,pttwo,
     & ptll(4),Etmiss(4),Etll,Rpt,pt,mtransZZ,mtransWW,mtransWZ
      integer:: switch,n,nplotmax,ilep,ijet,ielect,inu,
     & ide(4),idl(4),idj(4),idn(4),j
      integer tag
      logical:: is_lepton,is_hadronic,
     & is_neutrino
      logical, save::first=.true.
      common/nplotmax/nplotmax
      save ilep,ijet,inu,idl,idj,idn
ccccc!$omp threadprivate(first,/nplotmax/,ilep,ijet,inu,idl,idj,idn)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************
      
      if (first) then
        ilep=0
        ijet=0
        inu=0

C----setup particle identification 
        do j=3,8
          if (is_lepton(j)) then      ! lepton
            ilep=ilep+1
            idl(ilep)=j
          endif
          if (is_hadronic(j)) then      ! jet 
            ijet=ijet+1
            idj(ijet)=j
          endif
          if (is_neutrino(j)) then      ! neutrino
            inu=inu+1
            idn(inu)=j
          endif
        enddo
        
        if ((ilep < 2) .and. (removebr .eqv. .false.)) then
          write(6,*) 'Error in plotting routine:'
          write(6,*) 'did not find at least 2 leptons: ',ilep
          stop
        endif

        if (ijet < 2) then
          write(6,*) 'Error in plotting routine:'
          write(6,*) 'did not find at least 2 jets: ',ijet
          stop
        endif
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

c      if (ilep == 4) then 
      m3456=(p(idl(1),4)+p(idl(2),4)+p(idl(3),4)+p(idl(4),4))**2
     &     -(p(idl(1),1)+p(idl(2),1)+p(idl(3),1)+p(idl(4),1))**2
     &     -(p(idl(1),2)+p(idl(2),2)+p(idl(3),2)+p(idl(4),2))**2
     &     -(p(idl(1),3)+p(idl(2),3)+p(idl(3),3)+p(idl(4),3))**2
      m3456=sqrt(max(m3456,zip))
c      else 
c      m3456=0._dp
c      endif
c--- mtrans defined according to Eq.(13) in ATLAS-CONF-2014-042
      if (ilep == 2) then 
      ptll(:)=p(idl(1),:)+p(idl(2),:)
      Etmiss(:)=p(idn(1),:)+p(idn(2),:)
      mtransZZ=
     & (sqrt(max(zmass**2+ptll(1)**2+ptll(2)**2,zip))
     & +sqrt(max(zmass**2+Etmiss(1)**2+Etmiss(2)**2,zip)))**2
     &-(ptll(1)+Etmiss(1))**2-(ptll(2)+Etmiss(2))**2
      mtransZZ=sqrt(max(mtransZZ,zip))
      
      Etll=sqrt(max(ptll(4)**2-ptll(3)**2,zip))
      mtransWW=
     & (Etll+sqrt(Etmiss(1)**2+Etmiss(2)**2))**2
     &-(ptll(1)+Etmiss(1))**2-(ptll(2)+Etmiss(2))**2
      mtransWW=sqrt(max(mtransWW,zip))      
      
      pt34=pttwo(idl(1),idl(2),p)
      Rpt=pt(idl(1),p)*pt(idl(2),p)
     &  /(pt(idj(1),p)*pt(idj(2),p))
      else
      Rpt=zip
      mtransZZ=zip
      mtransWW=zip
      endif
c--- mtrans defined according to Eq.(6.5) in 1412.8367
      if (ilep == 3) then 
      ptll(:)=p(idl(1),:)+p(idl(2),:)+p(idl(3),:)
      Etmiss(:)=p(idn(1),:)
      mtransWZ=
     & (sqrt(max(ptll(4)**2-ptll(3)**2,zip))
     & +sqrt(max(Etmiss(1)**2+Etmiss(2)**2,zip)))**2
     &-(ptll(1)+Etmiss(1))**2-(ptll(2)+Etmiss(2))**2
      mtransWZ=sqrt(max(mtransWZ,zip))
      else
      mtransWZ=zip
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

c--- Plots of m(3456) in specific regions
      call bookplot(n,tag,'0 < m(3456) < 2000',
     & m3456,wt,wt2,zip,2000._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'100 < m(3456) < 600',
     & m3456,wt,wt2,100._dp,600._dp,5._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'600 < m(3456) < 1100',
     & m3456,wt,wt2,600._dp,1100._dp,5._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'1100 < m(3456) < 1600',
     & m3456,wt,wt2,1100._dp,1600._dp,5._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'0 < m(3456) < 100',
     & m3456,wt,wt2,zip,100._dp,5._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(3456) < 2010',
     & m3456,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < m(3456) < 2010',
     & m3456,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < m(3456) < 2020',
     & m3456,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(3456) < 130',
     & m3456,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'0 < mtransZZ < 2000',
     & mtransZZ,wt,wt2,zip,2000._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransZZ < 2010',
     & mtransZZ,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < mtransZZ < 2010',
     & mtransZZ,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < mtransZZ < 2020',
     & mtransZZ,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransZZ < 130',
     & mtransZZ,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'0 < mtransWW < 2000',
     & mtransWW,wt,wt2,0._dp,2000._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransWW < 2010',
     & mtransWW,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < mtransWW < 2010',
     & mtransWW,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < mtransWW < 2020',
     & mtransWW,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransWW < 130',
     & mtransWW,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'0 < mtransWZ < 2000',
     & mtransWZ,wt,wt2,0._dp,2000._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransWZ < 2010',
     & mtransWZ,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < mtransWZ < 2010',
     & mtransWZ,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < mtransWZ < 2020',
     & mtransWZ,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < mtransWZ < 130',
     & mtransWZ,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'Rpt',
     & Rpt,wt,wt2,0._dp,50._dp,0.5_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'pt(Z)',
     & pt34,wt,wt2,0._dp,2._dp,0.02_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'+INTEGRAL+ pt(Z)',
     & pt34,wt,wt2,0._dp,10._dp,0.1_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'50 < m(3456) < 250',
     & m3456,wt,wt2,50._dp,250._dp,2._dp,'log')
      n=n+1
      
c--- usual plots for 3+4
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for 5+6
      call autoplot2(p,56,5,6,tag,wt,wt2,n)

c--- usual plots for 3+4+5+6
      call autoplot4(p,3456,3,4,5,6,tag,wt,wt2,n)

c--- usual plots for 3+6
      call autoplot2(p,36,3,6,tag,wt,wt2,n)

c--- usual plots for 4+5
      call autoplot2(p,45,4,5,tag,wt,wt2,n)

c--- usual plots for 7+8
      call autoplot2(p,78,7,8,tag,wt,wt2,n)

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
