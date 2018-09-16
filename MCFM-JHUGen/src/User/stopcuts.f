      subroutine stopcuts(p,maxparts,ht,qeta,mlbnu,merecon,reconcorr)
c--- given the event momenta in p (with maxparts entries),
c--- performs the single-top search cuts and then returns the
c--- values of HT, Q*eta and M(l+b+nu) for the plotting routine
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'process.f'
      integer j,ilept,ibjet,ispect,maxparts,countjet,jetindex(mxpart)
      double precision pt,etarap,p(mxpart,4),ran2,
     . etvec(2),plept(4),etmissing(4),
     . mlbnu,ht,qeta,missinget,merecon,reconcorr
            
c--- initialize all variables      
      ilept=-1
      ibjet=-1
      ispect=-1
      ht=-1d0
      qeta=-99d0
      mlbnu=-1d0
      merecon=-1d0
      reconcorr=-99d0
      
c--- find one lepton that satisfies the pt and rapidity cuts
      do j=3,maxparts
        if (     (plabel(j).eq.'el') .or. (plabel(j).eq.'ea')
     .      .or. (plabel(j).eq.'ml') .or. (plabel(j).eq.'ma')) then
          if (       (pt(j,p) .gt. 20d0) .and.
     .      (abs(etarap(j,p)) .lt. 1.1d0)) then
            if (ilept .eq. -1) then
              ilept=j
            else
              goto 999
            endif
          endif
        endif
      enddo

c--- return if there are no leptons that satisfy the cuts      
      if (ilept .eq. -1) goto 999
      
c--- identify the jets
      countjet=0      
      do j=3,maxparts
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

      if (countjet .ne. 2) then
        write(6,*) 'Too many jets identified'
        stop
      endif

c--- the above is too inefficient for W+2 jets: apply factor after       
      if (case .eq. 'W_2jet') then
        if (pt(jetindex(1),p) .gt.  pt(jetindex(2),p)) then
c        if (r1 .gt. 0.5d0) then
          jetlabel(1)='bq'
          jetlabel(2)='pp'
        else
          jetlabel(1)='pp'
          jetlabel(2)='bq'
        endif
      endif             
             
c--- find the leading b-jet and the spectator jet    
      if ((jetlabel(1) .eq. 'bq') .or. (jetlabel(1) .eq. 'ba')) then
        ibjet=1
        ispect=2
      endif
      if ((jetlabel(2) .eq. 'bq') .or. (jetlabel(2) .eq. 'ba')) then
        if (ibjet .eq. -1) then
          ibjet=2
          ispect=1
        else    
c--- this is an arbitrary condition (leading jet?) for the case of 2 b's
c          if (pt(jetindex(1),p) .gt.  pt(jetindex(2),p)) then
          if (ran2() .lt. 0.5d0) then
            ibjet=1
            ispect=2
          else
            ibjet=2
            ispect=1
          endif  
        endif
      endif
       
      if (ibjet .eq. -1) goto 999 
       
c--- form the missing et vector 
      do j=1,4
        plept(j)=p(ilept,j)
        if ((j.eq.1) .or. (j.eq.2) )
     .  etvec(j)=-p(jetindex(ibjet),j)-p(jetindex(ispect),j)
     .           -p(ilept,j)
      enddo
      
      call wconstruct(etvec,plept,plabel(ilept),etmissing)

      missinget=dsqrt(etmissing(1)**2+etmissing(2)**2)
      
c--- Perform the cut on missing Et      
      if (missinget .lt. 20d0) goto 999
      
c--- variable to tell if reconstruction is correct or not
      if (abs(etmissing(4)-p(7-ilept,4)) .lt. 1d-6*etmissing(4)) then
        reconcorr=0.5d0
      else
        reconcorr=-0.5d0
      endif
      
      merecon=etmissing(4)
      
c--- form the invariant mass of the lepton, b-jet and constructed
c---  momentum of putative neutrino
      mlbnu=(p(ilept,4)+p(jetindex(ibjet),4)+etmissing(4))**2
      do j=1,3
        mlbnu=mlbnu
     .     -(p(ilept,j)+p(jetindex(ibjet),j)+etmissing(j))**2 
      enddo
      mlbnu=dsqrt(max(mlbnu,0d0))
c      write(6,*) ilept,jetindex(ibjet),inu
c      write(6,*) mlbnu
      
c--- Perform the cut on this invariant mass     
      if ((mlbnu .lt. 140d0) .or. (mlbnu .gt. 210d0)) then
        mlbnu=-1d0
        merecon=-1d0
        reconcorr=-99d0
        goto 999
      endif
      
      ht=pt(ilept,p)+missinget
     .  +pt(jetindex(1),p)+pt(jetindex(2),p)

c--- leading jet condition, only applied to qeta for t-channel search   
      if (max(pt(jetindex(1),p),pt(jetindex(2),p)) .gt. 30d0) then
        if ((plabel(ilept).eq.'el') .or. (plabel(ilept).eq.'ml')) then
          qeta=-1d0
c--- include e+ only for W^+
          mlbnu=-1d0
          merecon=-1d0
          reconcorr=-99d0
          ht=-1d0
          qeta=-99d0
          goto 999
        else
          qeta=+1d0
        endif
        qeta=qeta*etarap(jetindex(ispect),p)
      endif
      
      return
  
c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return    
      
      end
      
     
