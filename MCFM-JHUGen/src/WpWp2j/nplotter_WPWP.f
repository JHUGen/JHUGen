  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!
      subroutine nplotter_WpWp(p,wt,wt2,switch)
      implicit none
      include 'types.f'
      include 'vegas_common.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'outputflags.f'
      include 'nqcdjets.f'
      integer n,switch,i7,i8,i9,nu,nplotmax
      integer tag
      real(dp):: PT,WT,WT2,p(mxpart,4),fphi,etarap,etaraptwo
      real(dp):: tmp7(4),tmp8(4),tmp9(4)
      integer eventpart
      logical first,jetmerge
      common/nplotmax/nplotmax
      common/jetmerge/jetmerge
      data first/.true./
      save first,eventpart
      !------------
      real(dp):: tomphill, tommtww,tompt4,tompt5
      real(dp):: tometa4,tometa5,tometa45,tompt7
      real(dp):: tompt8, tometa7,tometa8,tometa78
      real(dp):: tompt9, tometa9, tometa79,tometa89
      real(dp):: tometmiss,tomHT,tomHTTOT,tomHTJET
      real(dp):: tomptmiss(4),tomptemu(4),tommemusq,tometemu
      real(dp):: rttommemusq,tomMetmiss
      save tomphill,tommtww,tompt4,tompt5,tometa4,tometa5,tometa45
      save tompt7,tometa7,tompt8,tometa8,tometa78,tompt9,tometa9
      save tometa79,tometa89,tomHTJET,tometmiss,tomHT,tomHTTOT
!$omp threadprivate(first,eventpart)
!$omp threadprivate(tomphill,tommtww,tompt4,tompt5,tometa4,tometa5)
!$omp threadprivate(tompt7,tometa7,tompt8,tometa8,tometa78,tompt9)
!$omp threadprivate(tometa79,tometa89,tomHTJET,tometmiss,tomHT,tomHTTOT)
!$omp threadprivate(tometa45,tometa9)
!$omp threadprivate(/jetmerge/)      
ccccc!$omp threadprivate(/nplotmax/)      

      if (first) then

        tag=tagbook

      !TOM added
        tomphill=zero
        tommtww=zero
        tompt4=zero
        tompt5=zero
        tometa4=zero
        tometa5=zero
        tometa45=zero
        tompt7=zero
        tometa7=zero
        tompt8=zero
        tometa8=zero
        tometa78=zero
        tompt9=zero
        tometa9=zero
        tometa79=zero
        tometa89=zero
        tomHTJET=zero
        tometmiss=zero
        tomHT=zero
        tomHTTOT=zero
        !END TOM ADDED

        jetmerge=.true.
        jets=nqcdjets
        eventpart = 4+jets
        goto 99
      else
        tag=tagplot
      endif

c--- 'eventpart' will contain the number of actual particles that have
c--- a defined momentum. For most processes, this is calculated as follows:
c----  for lowest order and virtual terms switch=0 and eventpart=npart+2
c---   for real events switch=0 and eventpart=npart+2
c---   for real counter-events switch=1 and eventpart=npart+1
c--- There are some processes for which this is not correct and these
c---  are handled with reference to nproc  
      eventpart=npart-switch+2

      if (jets .ge. 0) then
        eventpart=4+jets
      endif

         tompt7=zero
         tompt8=zero
         i7=7
         i8=8
         i9=9

         tompt7=pt(7,p)
         tompt8=pt(8,p)
         if (jets .eq. 3) tompt9=pt(9,p)

         !sort for 2jets
         if (jets .eq. 2) then
            if (tompt8.gt.tompt7) then
            i7=8
            i8=7
            endif
         endif
         !sort for 3jets
         if (jets .eq. 3) then
            if ((tompt7.gt.tompt8).and.(tompt7.gt.tompt9)) then
               i7=7
               if (tompt8.gt.tompt9) then
                  i8=8
                  i9=9
               else
                  i8=9
                  i9=8
               endif
            endif
            if ((tompt8.gt.tompt7).and.(tompt8.gt.tompt9)) then
               i7=8
               if (tompt7.gt.tompt9) then
                  i8=7
                  i9=9
               else
                  i8=9
                  i9=7
               endif
            endif
            if ((tompt9.gt.tompt7).and.(tompt9.gt.tompt8)) then
               i7=9
               if (tompt7.gt.tompt8) then
                  i8=7
                  i9=8
               else
                  i8=8
                  i9=7
               endif
            endif
         endif
           
         !---swap momentum if needed
         do nu=1,4
            tmp7(nu)=p(i7,nu)
            tmp8(nu)=p(i8,nu)
            tmp9(nu)=p(i9,nu)
         enddo
         do nu=1,4
            p(7,nu)=tmp7(nu)
            p(8,nu)=tmp8(nu)
            p(9,nu)=tmp9(nu)
         enddo

         !Setup interesting variables
         !antilepton
         tompt4=pt(4,p)
         tometa4=etarap(4,p)
         !lepton
         tompt5=pt(6,p)
         tometa5=etarap(6,p)
         !rapidity and angle between leptons
         tometa45=etaraptwo(4,6,p)
         tomphill=fphi(4,6,p)
         !missing energy
         do nu=1,4
         tomptmiss(nu) = p(3,nu)+p(5,nu)
         tomptemu(nu) = p(4,nu)+p(6,nu)
         enddo


         !jet rapidities
         tompt7=pt(7,p)
         tometa7=etarap(7,p)
         tompt8=pt(8,p)
         tometa8=etarap(8,p)
         tometa78=etaraptwo(7,8,p)


         if (jets .eq. 3) then
         tompt9=pt(9,p)
         tometa9=etarap(9,p)
         tometa79=etaraptwo(7,9,p)
         tometa89=etaraptwo(8,9,p)
         endif

c         jets210 = log_val_opt('-jets210',.false.)
c
c         if (jets210) then
c
c         if (jets.ge.2) then
c         if ((tompt7 .gt. 30._dp).and.(tompt8.gt.30._dp)) jets210int30=2
c         if ((tompt7 .gt. 30._dp).and.(tompt8.lt.30._dp)) jets210int30=1
c         if ((tompt7 .lt. 30._dp).and.(tompt8.lt.30._dp)) jets210int30=0
c         endif
c         if (jets.eq.3) then
c            if (tompt9.gt.30._dp) jets210int30=3
c         endif
c
c         if (jets210int30.ge.2) then
c            tomRlj1= R(p,4,7)
c            tomRlj2= R(p,4,8)
c         endif
c         endif

         tometmiss=sqrt( tomptmiss(1)**2+tomptmiss(2)**2)

         !transverse WW mass
         tommemusq=tomptemu(4)**2-tomptemu(1)**2-tomptemu(2)**2
     .                -tomptemu(3)**2
         tomMetmiss=sqrt( tomptmiss(1)**2+tomptmiss(2)**2 + tommemusq )
         tometemu = sqrt( tomptemu(1)**2+tomptemu(2)**2+tommemusq )
         tommtww=sqrt( (tomMetmiss+tometemu)**2 - 
     .                     (tomptmiss(1)+tomptemu(1))**2 -
     .                     (tomptmiss(2)+tomptemu(2))**2 )

         !HT 
         tomHTJET = tompt7+tompt8
         tomHT= tompt7+tompt8+tometmiss
         tomHTTOT = tompt4+tompt5+tompt7+tompt8+tometmiss
         
         if (jets .eq. 3) then
         tomHTJET = tomHTJET+tompt9
         tomHT = tomHT + tompt9
         tomHTTOT = tomHTTOT + tompt9
         endif

  99  continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/real(itmx,dp))  
      endif

      n=1                        

      !these histograms are for no jet210 cuts appied to jets
      call bookplot(n,tag,'HTTOT',tomHTTOT,
     .                      wt,wt2,zero,4000._dp,50._dp,'lin')
      n=n+1      
      call bookplot(n,tag,'HT',tomHT,
     .                      wt,wt2,zero,3000._dp,30._dp,'lin')
      n=n+1            
      call bookplot(n,tag,'HTJET',tomHTJET,
     .                      wt,wt2,zero,3000._dp,30._dp,'lin')
      n=n+1            
      call bookplot(n,tag,'pt_l1',tompt4,wt,wt2,zero,450._dp,5._dp,'log')
      n=n+1      
      call bookplot(n,tag,'pt_l2',tompt5,wt,wt2,zero,450._dp,5._dp,'log')
      n=n+1      
      call bookplot(n,tag,'eta_l1',tometa4,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      
      call bookplot(n,tag,'eta_l2',tometa5,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      
      call bookplot(n,tag,'eta_ll',tometa45,
     .                       wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      
      call bookplot(n,tag,'pt_j1',tompt7,wt,wt2,zero,800._dp,10._dp,'log')
      n=n+1      
      call bookplot(n,tag,'pt_j2',tompt8,wt,wt2,zero,450._dp,5._dp,'log')
      n=n+1      
      call bookplot(n,tag,'eta_j1',tometa7,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      
      call bookplot(n,tag,'eta_j2',tometa8,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      
      call bookplot(n,tag,'eta_j1j2',tometa78,
     .           wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1      

      if (eventpart .gt. 6) then
      call bookplot(n,tag,'pt_j3',tompt9,wt,wt2,zero,270._dp,3._dp,'log')
      endif
      n=n+1      
      if (eventpart .gt. 6) then
      call bookplot(n,tag,'eta_j3',tometa9,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      endif
      n=n+1
      if (eventpart .gt. 6) then
      call bookplot(n,tag,'eta_j1j3',tometa79,
     .           wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      endif
      n=n+1      
      if (eventpart .gt. 6) then
      call bookplot(n,tag,'eta_j2j3',tometa89,
     .           wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      endif
      n=n+1
      call bookplot(n,tag,'Etmiss',tometmiss,
     .            wt,wt2,zero,800._dp,10._dp,'lin')
      n=n+1         
      call bookplot(n,tag,'mtww',tommtww,
     .                   wt,wt2,zero,2500._dp,50._dp,'lin')
      n=n+1      
      rttommemusq=sqrt(tommemusq)
      call bookplot(n,tag,'mll',rttommemusq,
     .                   wt,wt2,zero,1200._dp,30._dp,'lin')
      n=n+1      
      call bookplot(n,tag,'Phi_ll',tomphill,
     .                      wt,wt2,zero,3.5_dp,0.1_dp,'lin')
      n=n+1      

      n=n-1

      if (n .gt. maxhisto) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > ',maxhisto,', which is the built-in maximum.'
        write(6,*) 'To use more histograms, change the value of the'
        write(6,*) 'constant MAXHISTO in src/Inc/nplot.f then do:'
        write(6,*)
        write(6,*) ' make clean; make        to recompile from scratch.'
        write(6,*)
        stop
      endif
c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return 
      end
