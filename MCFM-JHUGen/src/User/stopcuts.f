      subroutine stopcuts(p,maxparts,ht,qeta,mlbnu,merecon,reconcorr)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with maxparts entries),
c--- performs the single-top search cuts and then returns the
c--- values of HT, Q*eta and M(l+b+nu) for the plotting routine
      include 'constants.f'
      include 'mxpart.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'kprocess.f'
      integer:: j,ilept,ibjet,ispect,maxparts,countjet,jetindex(mxpart)
      real(dp):: pt,etarap,p(mxpart,4),ran2,
     & etvec(2),plept(4),etmissing(4),
     & mlbnu,ht,qeta,missinget,merecon,reconcorr
            
c--- initialize all variables      
      ilept=-1
      ibjet=-1
      ispect=-1
      ht=-1._dp
      qeta=-99._dp
      mlbnu=-1._dp
      merecon=-1._dp
      reconcorr=-99._dp
      
c--- find one lepton that satisfies the pt and rapidity cuts
      do j=3,maxparts
        if (     (plabel(j)=='el') .or. (plabel(j)=='ea')
     &      .or. (plabel(j)=='ml') .or. (plabel(j)=='ma')) then
          if (       (pt(j,p) > 20._dp) .and.
     &      (abs(etarap(j,p)) < 1.1_dp)) then
            if (ilept == -1) then
              ilept=j
            else
              goto 999
            endif
          endif
        endif
      enddo

c--- return if there are no leptons that satisfy the cuts      
      if (ilept == -1) goto 999
      
c--- identify the jets
      countjet=0      
      do j=3,maxparts
        if (     (plabel(j) == 'pp') .or. (plabel(j) == 'qj')
     &      .or. (plabel(j) == 'bq') .or. (plabel(j) == 'ba')) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

      if (countjet .ne. 2) then
        write(6,*) 'Too many jets identified'
        stop
      endif

c--- the above is too inefficient for W+2 jets: apply factor after       
      if (kcase==kW_2jet) then
        if (pt(jetindex(1),p) >  pt(jetindex(2),p)) then
c        if (r1 > half) then
          jetlabel(1)='bq'
          jetlabel(2)='pp'
        else
          jetlabel(1)='pp'
          jetlabel(2)='bq'
        endif
      endif             
             
c--- find the leading b-jet and the spectator jet    
      if ((jetlabel(1) == 'bq') .or. (jetlabel(1) == 'ba')) then
        ibjet=1
        ispect=2
      endif
      if ((jetlabel(2) == 'bq') .or. (jetlabel(2) == 'ba')) then
        if (ibjet == -1) then
          ibjet=2
          ispect=1
        else    
c--- this is an arbitrary condition (leading jet?) for the case of 2 b's
c          if (pt(jetindex(1),p) >  pt(jetindex(2),p)) then
          if (ran2() < half) then
            ibjet=1
            ispect=2
          else
            ibjet=2
            ispect=1
          endif  
        endif
      endif
       
      if (ibjet == -1) goto 999 
       
c--- form the missing et vector 
      do j=1,4
        plept(j)=p(ilept,j)
        if ((j==1) .or. (j==2) )
     &  etvec(j)=-p(jetindex(ibjet),j)-p(jetindex(ispect),j)
     &           -p(ilept,j)
      enddo
      
      call wconstruct(etvec,plept,plabel(ilept),etmissing)

      missinget=sqrt(etmissing(1)**2+etmissing(2)**2)
      
c--- Perform the cut on missing Et      
      if (missinget < 20._dp) goto 999
      
c--- variable to tell if reconstruction is correct or not
      if (abs(etmissing(4)-p(7-ilept,4)) < 1.e-6_dp*etmissing(4)) then
        reconcorr=half
      else
        reconcorr=-half
      endif
      
      merecon=etmissing(4)
      
c--- form the invariant mass of the lepton, b-jet and constructed
c---  momentum of putative neutrino
      mlbnu=(p(ilept,4)+p(jetindex(ibjet),4)+etmissing(4))**2
      do j=1,3
        mlbnu=mlbnu
     &     -(p(ilept,j)+p(jetindex(ibjet),j)+etmissing(j))**2 
      enddo
      mlbnu=sqrt(max(mlbnu,zip))
c      write(6,*) ilept,jetindex(ibjet),inu
c      write(6,*) mlbnu
      
c--- Perform the cut on this invariant mass     
      if ((mlbnu < 140._dp) .or. (mlbnu > 210._dp)) then
        mlbnu=-1._dp
        merecon=-1._dp
        reconcorr=-99._dp
        goto 999
      endif
      
      ht=pt(ilept,p)+missinget
     &  +pt(jetindex(1),p)+pt(jetindex(2),p)

c--- leading jet condition, only applied to qeta for t-channel search   
      if (max(pt(jetindex(1),p),pt(jetindex(2),p)) > 30._dp) then
        if ((plabel(ilept)=='el') .or. (plabel(ilept)=='ml')) then
          qeta=-1._dp
c--- include e+ only for W^+
          mlbnu=-1._dp
          merecon=-1._dp
          reconcorr=-99._dp
          ht=-1._dp
          qeta=-99._dp
          goto 999
        else
          qeta=+1._dp
        endif
        qeta=qeta*etarap(jetindex(ispect),p)
      endif
      
      return
  
c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return    
      
      end
      
     
