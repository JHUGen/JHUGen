      subroutine hwwcuts(p,maxparts,dphi_ll,m_ll,mtrans,scut1,scut2)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic H->WW search cuts and then returns
c--- the values of dphi_ll,m_ll,mtrans,scut1,scut2 for the
c--- plotting routine; a point that fails the cuts returns mtrans=-1d0
      implicit none
      include 'constants.f'
      include 'plabel.f'
      integer ilept(2),ic,j,maxparts
      double precision pt,etarap,p(mxpart,4),etvec(2),plept(4),
     . missinget,dphi_ll,m_ll,mtrans,scut1,scut2,cosdphi,maxpt,minpt
            
c--- initialize all variables  
      dphi_ll=-1d0
      m_ll=-1d0
      mtrans=-1d0
      scut1=-1d0
      scut2=-1d0    
      etvec(1)=0d0
      etvec(2)=0d0

************************** START BASIC CUTS ****************************

      ic=0      
c--- find two leptons that satisfy the pt and rapidity cuts
      do j=3,maxparts
        if (     (plabel(j).eq.'el') .or. (plabel(j).eq.'ea')
     .      .or. (plabel(j).eq.'ml') .or. (plabel(j).eq.'ma')) then
          if (       (pt(j,p) .gt. 20d0) .and.
     .      (abs(etarap(j,p)) .lt. 2.5d0)) then
            ic=ic+1
            ilept(ic)=j
          endif
        endif
      enddo

c--- return if there are not 2 leptons that satisfy the cuts      
      if (ic .ne. 2) goto 999      
      
c--- form the combined lepton momentum and the missing et vector 
      do j=1,4
        plept(j)=p(ilept(1),j)+p(ilept(2),j)
        if ((j.eq.1) .or. (j.eq.2))
     .  etvec(j)=etvec(j)-p(ilept(1),j)-p(ilept(2),j)
      enddo

      missinget=dsqrt(etvec(1)**2+etvec(2)**2)
      
c--- perform the cut on missing Et, rejecting less than 30 GeV   
      if (missinget .lt. 25d0) goto 999
     
c---  perform veto if the event contains jets of 30 GeV with |eta|<3
      do j=3,maxparts
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          if (       (pt(j,p) .gt. 20d0) .and.
     .       (abs(etarap(j,p)) .lt. 3d0)) then
            goto 999
          endif
        endif
      enddo
      
c--- opening angle between the leptons in the transverse plane
      dphi_ll=
     .   (p(ilept(1),1)*p(ilept(2),1)+p(ilept(1),2)*p(ilept(2),2))
     .   /dsqrt((p(ilept(1),1)**2+p(ilept(1),2)**2)
     .         *(p(ilept(2),1)**2+p(ilept(2),2)**2))
      if (dphi_ll .lt. -0.999999999D0) dphi_ll=-1d0
      dphi_ll=dacos(dphi_ll) 
      
      if (dphi_ll .gt. pi/4d0) then
        goto 999
      endif
            
c--- dilepton invariant mass
      m_ll=dsqrt(plept(4)**2-plept(1)**2-plept(2)**2-plept(3)**2)

      if (m_ll .gt. 35d0) then
        goto 999
      endif
      
      maxpt=max(pt(ilept(1),p),pt(ilept(2),p))
      minpt=min(pt(ilept(1),p),pt(ilept(2),p))
            
      if (minpt .lt. 25d0) then
        goto 999
      endif
      
      if ((maxpt .lt. 35d0) .or. (maxpt .gt. 50d0)) then
        goto 999
      endif
      
*************************** END BASIC CUTS *****************************

******************** CALCULATE QUANTITIES TO PLOT **********************

c--- opening angle between the leptons in the transverse plane
      dphi_ll=
     .   (p(ilept(1),1)*p(ilept(2),1)+p(ilept(1),2)*p(ilept(2),2))
     .   /dsqrt((p(ilept(1),1)**2+p(ilept(1),2)**2)
     .         *(p(ilept(2),1)**2+p(ilept(2),2)**2))
      if (dphi_ll .lt. -0.999999999D0) dphi_ll=-1d0
      dphi_ll=dacos(dphi_ll) 
            
c--- dilepton invariant mass
      m_ll=dsqrt(plept(4)**2-plept(1)**2-plept(2)**2-plept(3)**2)

c--- transverse mass
c---    first compute the cosine of the azimuthal angle between the
c---    lepton system and the missing Et (often, this will simply be pi)
      cosdphi=
     .   (plept(1)*etvec(1)+plept(2)*etvec(2))
     .   /dsqrt((plept(1)**2+plept(2)**2)
     .         *(etvec(1)**2+etvec(2)**2))
c---    transverse mass calculation
      mtrans=2d0*dsqrt(plept(1)**2+plept(2)**2)*missinget*(1d0-cosdphi)
      mtrans=dsqrt(max(mtrans,0d0))

********************** END OF QUANTITIES TO PLOT ***********************

********************* START OF CUT CROSS SECTIONS **********************

c--- values of scut1,scut2 less than zero indicate that this point
c--- should not contribute to those cross sections

c--- require the azimuthal angle to be less than 45 degrees
      if (dphi_ll .gt. pi/4d0) goto 999   
      
c--- require the dilepton mass to be less than 35 GeV       
      if (m_ll .gt. 35d0) goto 999                  
      
c--- scut1 corresponds to the cross section with 125 < mtrans < 155
      if ((mtrans .ge. 125d0) .and. (mtrans .le. 155d0)) scut1=1d0    
      
c--- scut2 corresponds to the cross section with 140 < mtrans < 180
      if ((mtrans .ge. 140d0) .and. (mtrans .le. 180d0)) scut2=1d0    
      
********************** END OF CUT CROSS SECTIONS ***********************

      return
  

c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return    
      
      end
      
     
