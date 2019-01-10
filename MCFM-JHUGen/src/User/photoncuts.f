      logical function photoncuts(pjet)
************************************************************************
*   Author: J.M. Campbell, 24th January 2011                           *
*       and C. Williams                                                *
*   This routine imposes the photon cuts that are specified in the     *
*   input file. The cuts are applied here (rather than in gencuts.f)   *
*   since they may be necessary for a finite cross section             *
*   and thus should be applied regardless of the value of "makecuts"   *
*                                                                      *
*   Return TRUE if this point FAILS the cuts                           *
*   Modified March 11 to apply mass cuts to m34 for gamgam (CW)        *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'process.f'
      include 'leptcuts.f'
      include 'z_dip.f'
      include 'jetlabel.f'
      logical first,is_lepton,is_photon,is_hadronic
      integer i,j,k,im1,im2
      double precision ptm,ptm1,ptm2,pt3,pt4
      integer countlept,leptindex(mxpart),countgamm,gammindex(mxpart),
     & countjet,jetindex(mxpart)
      double precision pjet(mxpart,4),pt,etarap,R,pt1,pt2,pth,pts,s34
      data first/.true./
      save first,countlept,countgamm,countjet,
     & leptindex,gammindex,jetindex
!$omp threadprivate(first,countlept,countgamm,countjet)
!$omp threadprivate(leptindex,gammindex,jetindex)
      
      photoncuts=.false.

      if (first) then
      first=.false.      
c--- write-out the cuts we are using
      write(6,*)
      write(6,*)  '****************** Photon cuts *********************'
      write(6,*)  '*                                                  *'
      write(6,99) '*   pt(photon 1)         >   ',gammpt,
     &                '                *'
      write(6,99) '*   pt(photon 2)         >   ',gammpt2,
     &                '                *'
      if((case.eq.'trigam') .or. (case.eq.'fourga'))then 
         write(6,99) '*   pt(photon 3)         >   ',gammpt3,
     &        '                *'
      endif
      if(case.eq.'fourga') then 
         write(6,99) '*   pt(photon 4)         >   ',gammpt3,
     &        '                *'
      endif
      write(6,99) '*   eta(photon)          <   ',gammrap,
     &                '                *'
      write(6,99) '*   R(photon,lepton)     >   ',Rgalmin,
     &                '                *'
      write(6,99) '*   R(photon,photon)     >   ',Rgagamin,
     &                '                *'
      write(6,99) '*   R(photon,jet)        >   ',Rgajetmin,
     &                '                *'
      write(6,*)  '*                                                  *'
      write(6,*)  '****************************************************'
      if(case.eq.'gamgam') then 
       write(6,*)
       write(6,*) '************* M(gam,gam) mass cuts *****************'
       write(6,*) '*                                                  *'
       write(6,98) dsqrt(wsqmin),'m34',dsqrt(wsqmax)
       write(6,*) '****************************************************'
      endif
c--- initialize counters and arrays that will be used to perform cuts      
      countlept=0
      countgamm=0
      countjet=0
      do j=3,mxpart
        if (is_lepton(j)) then
          countlept=countlept+1
          leptindex(countlept)=j
        endif
        if (is_photon(j)) then
          countgamm=countgamm+1
          gammindex(countgamm)=j
        endif
        if (is_hadronic(j)) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo
      endif      
      
C     Basic pt and rapidity cuts for photon
      if (countgamm .eq. 1) then
          if (     (pt(gammindex(1),pjet) .lt. gammpt) .or.
     &    (abs(etarap(gammindex(1),pjet)) .gt. gammrap)) then
            photoncuts=.true.
            return
          endif
      endif
      if (countgamm .eq. 2) then
        pt1=pt(gammindex(1),pjet) 
        pt2=pt(gammindex(2),pjet) 
        pth=max(pt1,pt2)
        pts=min(pt1,pt2)
        if ( ( pth .lt. gammpt) .or.
     &       ( pts .lt. gammpt2) .or.
     &       (abs(etarap(gammindex(1),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(2),pjet)) .gt. gammrap) ) then
          photoncuts=.true.
          return
        endif
      endif

      if (countgamm .eq. 3) then
        pt1=pt(gammindex(1),pjet) 
        pt2=pt(gammindex(2),pjet)
        pt3=pt(gammindex(3),pjet)
        pth=max(pt1,pt2,pt3)
        pts=min(pt1,pt2,pt3)
        ptm=pt1+pt2+pt3-pts-pth
        if ( ( pth .lt. gammpt) .or.
     &       ( ptm .lt. gammpt2) .or.
     &       ( pts .lt. gammpt3) .or.
     &       (abs(etarap(gammindex(1),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(2),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(3),pjet)) .gt. gammrap) ) then
          photoncuts=.true.
          return
        endif
      endif

      if (countgamm .eq. 4) then
        pt1=pt(gammindex(1),pjet) 
        pt2=pt(gammindex(2),pjet)
        pt3=pt(gammindex(3),pjet)
        pt4=pt(gammindex(4),pjet)
        pth=max(pt1,pt2,pt3,pt4)
        pts=min(pt1,pt2,pt3,pt4)
        do i=1,4 
           if((pt(gammindex(i),pjet).ne.pth)
     &          .and.(pt(gammindex(i),pjet).ne.pts)) then
           im1=i
           endif
        enddo
        do i=1,4 
           if((pt(gammindex(i),pjet).ne.pth)
     &          .and.(pt(gammindex(i),pjet).ne.pts).and.(i.ne.im1)) then
           im2=i
        endif
        enddo
        ptm1=max(pt(gammindex(im1),pjet),pt(gammindex(im2),pjet))
        ptm2=min(pt(gammindex(im1),pjet),pt(gammindex(im2),pjet))
        if ( ( pth .lt. gammpt) .or.
     &       ( ptm1 .lt. gammpt2) .or.
     &       ( ptm2 .lt. gammpt3) .or.
     &       ( pts .lt. gammpt3) .or.
     &       (abs(etarap(gammindex(1),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(2),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(3),pjet)) .gt. gammrap) .or.
     &       (abs(etarap(gammindex(4),pjet)) .gt. gammrap) ) then
          photoncuts=.true.
          return
        endif
      endif


c--- lepton-photon separation 
      if ((countlept .ge. 1) .and. (countgamm .ge. 1)) then
        do j=1,countgamm
        do k=1,countlept
          if (R(pjet,gammindex(j),leptindex(k)) .lt. Rgalmin) then
            photoncuts=.true.
            return
          endif
        enddo
        enddo
      endif
      
c--- photon-photon separation 
      if (countgamm .ge. 2) then
        do j=1,countgamm
        do k=j+1,countgamm
          if (R(pjet,gammindex(j),gammindex(k)) .lt. Rgagamin) then
            photoncuts=.true.
            return
          endif
        enddo
        enddo
      endif

c--- jet-photon separation (if there are 1 or more jets and photons)
      if ((jets .ge. 1) .and. (countgamm .ge. 1)) then
        do j=1,countgamm
        do k=1,jets
          if (R(pjet,gammindex(j),jetindex(k)) .lt. Rgajetmin) then
            photoncuts=.true.
            return
          endif
        enddo
        enddo
      endif

!--- Apply mass cuts here (gamgam only)
      if(case.eq.'gamgam') then 
         s34=+(pjet(3,4)+pjet(4,4))**2-(pjet(3,1)+pjet(4,1))**2
     &       -(pjet(3,2)+pjet(4,2))**2-(pjet(3,3)+pjet(4,3))**2
         if ((s34 .lt. wsqmin) .or. (s34 .gt. wsqmax)) then         
            photoncuts=.true. 
            return          
         endif
      endif
     



      return




c--- Lines below here are commented out for now;
c---  may be restored at a later date.
       
       
c--- jet-photon separation (if there are 1 or more jets and photons)
c      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
c        do j=1,countgamm
c        do k=1,njets
c          if (R(pjet,gammindex(j),jetindex(k)) .lt. Rjlmin) then
c            gencuts=.true.
c            return
c          endif
c        enddo
c        enddo
c      endif

c--- DEBUG: removed all isolation      
cc--- photon/hadron isolation     
c      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
c        do j=1,countgamm
cc--- Frixione cut, hep-ph/9801442
c          if (njets .gt. 2) then
c            write(6,*) 'Photon-hadron isolation not coded for njets > 2'
c            stop
c          endif
c          do k=1,njets
c            delta(k)=R(pjet,gammindex(j),jetindex(k))
c            ptjet(k)=pt(jetindex(k),pjet)
c          enddo
c          pntr=1
c          if (delta(2) .lt. delta(1)) pntr=2
c          if (njets .eq. 1) delta(2)=gammcone   
c          discr=(1d0-dcos(delta(pntr)))/(1d0-dcos(gammcone))
c          if (ptjet(pntr) .gt. discr*pt(gammindex(j),pjet)) then
c            gencuts=.true.
c            return
c          endif 
c          if (njets .ge. 2) then
c            discr=(1d0-dcos(delta(3-pntr)))/(1d0-dcos(gammcone))
c            if (ptjet(1)+ptjet(2) .gt. discr*pt(gammindex(j),pjet))then
c              gencuts=.true.
c              return
c            endif
c          endif
          
c--- this block was already removed
c--- optional jet-veto
c          if (   (pt(4+countgamm+1,pjet) .gt. 50d0)    
c     .     .and. (abs(etarap(4+countgamm+1,pjet)) .lt. 2.5d0)) then
c            gencuts=.true.
c          endif                     
c--- de-Florian,Signer cut
c          do nu=1,2
c            sumjetpt(nu)=0d0
c          enddo
c          do k=1,njets
c            if (R(pjet,gammindex(j),4+countgamm+k) .lt. gammcone) then
c              do nu=1,2
c              sumjetpt(nu)=sumjetpt(nu)+pjet(4+countgamm+k,nu)
c              enddo
c            endif
c          enddo
c          if ( dsqrt(sumjetpt(1)**2+sumjetpt(2)**2) 
c     .    .gt. gammcut*pt(gammindex(j),pjet) ) then
c            gencuts=.true.
c          endif
c--- this block was already removed

c        enddo
c      endif  
c--- DEBUG: removed all isolation      
      
 99   format(1x,a29,f6.2,a17)
 98   format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')
      end
 
 
 
 
