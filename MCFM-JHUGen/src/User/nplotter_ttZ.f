      subroutine nplotter_ttZ(p,wt,wt2,switch)
      implicit none
      include 'types.f'
c -- Taylored to PROCESS 532:
c -- pp -> t(-->nu(p3)+mu^+(p4)+b(p5))+t~(-->q(p7)+q~(p8)+b~(p6))+Z(e(p9),e~(p10))
c -- i.e. t decays semi-leptonically, tbar decays hadronically (nproc=532)
c -- Should be modified if t decays hadronically and tbar semi-leptonically (nproc=533)
c -- Assumes different flavor of the lepton from top decay to those from Z decay
c -- R. Rontsch 2013-03-01
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
      include 'jetlabel.f'
      include 'outputflags.f'
      real(dp):: p(mxpart,4),wt,wt2,tiny
      integer:: switch,n,nplotmax,jetindex(mxpart),ijet,
     & j
      real(dp):: pt,etarap, twomass,pttwo,
     &     ptthree,etarapthree,R
      real(dp):: wtjet(6),wtjet2(6)
      real(dp):: ptmiss,ptmu,ptem,ptep,ptj1,ptj2,ptj3,ptj4,
     & etamu,etaem,etaep,etaj1,etaj2,etaj3,etaj4,mZ,ptZ,
     & ptt,pttbar,etat,etatbar,DRttbar,DRtZ,DRtbarZ
      real(dp):: pb4dk(mxpart,4)
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

c -- jet-binned weights
       wtjet(:)=0._dp
       wtjet2(:)=0._dp
      if (jets == 0) then
         wtjet(1)=wt
         wtjet2(1)=wt2
      elseif (jets == 1) then
         wtjet(2)=wt
         wtjet2(2)=wt2
      elseif (jets == 2) then
         wtjet(3)=wt
         wtjet2(3)=wt2
      elseif (jets == 3) then
         wtjet(4)=wt
         wtjet2(4)=wt2
      elseif (jets == 4) then
         wtjet(5)=wt
         wtjet2(5)=wt2
      elseif (jets == 5) then
         wtjet(6)=wt
         wtjet2(6)=wt2
      endif

c -- identify the jets
      ijet=0
       do j=3,11
         if ((plabel(j) == 'pp') .or. (plabel(j) == 'bq')
     &   .or.(plabel(j) == 'ba')) then
           if (p(j,4) > tiny) then
             ijet=ijet+1
             jetindex(ijet)=j
           endif
         endif
       enddo   

c -- transverse momenta
      ptmiss=pt(3,p)
      ptmu=pt(4,p) 
      ptem=pt(9,p) 
      ptep=pt(10,p) 
      ptj1=pt(jetindex(1),p)
      ptj2=pt(jetindex(2),p)
      ptj3=pt(jetindex(3),p)
      ptj4=pt(jetindex(4),p)
      
c -- pseudorapidities
      etamu=etarap(4,p)
      etaem=etarap(9,p)
      etaep=etarap(10,p)
      etaj1=etarap(jetindex(1),p)
      etaj2=etarap(jetindex(2),p)
      etaj3=etarap(jetindex(3),p)
      etaj4=etarap(jetindex(4),p)

c -- quantities for Z:
      mZ=twomass(9,10,p)
      ptZ=pttwo(9,10,p)

c -- naive top properties
      ptt=ptthree(3,4,5,p)
      pttbar=ptthree(6,7,8,p)
      etat=etarapthree(3,4,5,p)
      etatbar=etarapthree(6,7,8,p)
c -- construct system prior to decay of tops and Z
      pb4dk(1,:)=p(3,:)+p(4,:)+p(5,:)
      pb4dk(2,:)=p(6,:)+p(7,:)+p(8,:)
      pb4dk(3,:)=p(9,:)+p(10,:)

c -- Delta R between tops and Z
      DRttbar=R(pb4dk,1,2)
      DRtZ=R(pb4dk,1,3)
      DRtbarZ=R(pb4dk,2,3)
      
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

c -- plot jet-binned cross-sections
      call bookplot(n,tag,'0jet',0.5_dp,wtjet(1),wtjet2(1),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'1jet',0.5_dp,wtjet(2),wtjet2(2),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'2jet',0.5_dp,wtjet(3),wtjet2(3),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'3jet',0.5_dp,wtjet(4),wtjet2(4),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'4jet',0.5_dp,wtjet(5),wtjet2(5),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'5jet',0.5_dp,wtjet(6),wtjet2(6),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1

c -- kinematic distributions
      call bookplot(n,tag,'pt,miss',ptmiss,
     &     wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(mu)',ptmu,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(e-)',ptem,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(e+)',ptep,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(j1)',ptj1,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(j2)',ptj2,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(j3)',ptj3,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(j4)',ptj4,wt,wt2,0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(mu)',etamu,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(e-)',etaem,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(e+)',etaep,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(j1)',etaj1,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(j2)',etaj2,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
       call bookplot(n,tag,'eta(j3)',etaj3,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
       call bookplot(n,tag,'eta(j4)',etaj4,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m(Z)',mZ,wt,wt2,0._dp,300._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(Z)',ptZ,wt,wt2,0._dp,300._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(t)',ptt,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt(tbar)',pttbar,
     &     wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(t)',etat,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(tbar)',etatbar,wt,wt2,
     &     -3.0_dp,3.0_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(t,tbar)',DRttbar,
     & wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1   
      call bookplot(n,tag,'DR(t,Z)',DRtZ,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1   
      call bookplot(n,tag,'DR(tbar,Z)',DRtbarZ,
     &     wt,wt2,0._dp,5._dp,0.1_dp,'lin')
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
