      subroutine cdfhwwcuts(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs H->WW search cuts along the lines of CDF Note 9887;
c--- a point that fails the cuts returns passed=.false.
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'plabel.f'
      integer ilept(2),ile,ineut(2),inu,i,j,maxparts
      logical passed
      double precision pt,etarap,p(mxpart,4),etvec(2),
     . missinget,m_ll,R,pt1,pt2,eta1,eta2,phi,phimin
c--- pt and eta thresholds for trigger lepton
      double precision trigpt,trigeta
c--- pt and eta thresholds for second lepton (should be looser)
      double precision scndpt,scndeta
c--- minimum threshold for dilepton invariant mass
      double precision m_llcut
c--- minimum threshold for missing Et "star" or "spec" (see defn below)
      double precision missingetcut
c--- size of cone for lepton isolation
      double precision Risol
c--- fractional pt threshold for lepton isolation
      double precision threshisol
      double precision mtmax,mtmin,mt45,m45,et_vec(4),etmiss,pttwo
      data trigpt,trigeta/20d0,0.8d0/
      data scndpt,scndeta/10d0,1.1d0/
      data m_llcut/16d0/
      data missingetcut/25d0/
      data Risol/0.4d0/
      data threshisol/0.1d0/

c--- default behaviour is that the cuts have been passed
      passed=.true.
                
c--- initialize all variables  
      etvec(1)=0d0
      etvec(2)=0d0

************************** START BASIC CUTS ****************************

      ile=0
      inu=0      
c--- identify the two leptons and neutrinos
      do j=3,maxparts
        if (     (plabel(j).eq.'el') .or. (plabel(j).eq.'ea')
     .      .or. (plabel(j).eq.'ml') .or. (plabel(j).eq.'ma')) then
          ile=ile+1
          ilept(ile)=j
        elseif  ((plabel(j).eq.'nl') .or. (plabel(j).eq.'na')) then
          inu=inu+1
          ineut(inu)=j
        endif
      enddo

c--- print warning if there are not 2 leptons     
      if (ile .ne. 2) then
        write(6,*) 'WARNING: cdfhwwcuts.f must be upgraded, <2 leptons'
      stop
      endif     
      
      
c--- compute lepton pt and rapidity
      pt1=pt(ilept(1),p)      
      pt2=pt(ilept(2),p)      
      eta1=abs(etarap(ilept(1),p))     
      eta2=abs(etarap(ilept(2),p))    

c--- CUT on basic acceptance cuts for the leptons
      if ((  (pt1 .gt. trigpt) .and. (eta1 .lt. trigeta)
     .  .and.(pt2 .gt. scndpt) .and. (eta2 .lt. scndeta)) .or.
     .    (  (pt2 .gt. trigpt) .and. (eta2 .lt. trigeta)
     .  .and.(pt1 .gt. scndpt) .and. (eta1 .lt. scndeta)) ) then
        continue ! cut passed
      else
        passed=.false.
      return
      endif

c--- compute invariant mass of leptons
      m_ll=dsqrt((p(ilept(1),4)+p(ilept(2),4))**2
     .          -(p(ilept(1),1)+p(ilept(2),1))**2
     .          -(p(ilept(1),2)+p(ilept(2),2))**2
     .          -(p(ilept(1),3)+p(ilept(2),3))**2)

c--- CUT on dilepton invariant mass
      if (m_ll .lt. m_llcut) then
        passed=.false.
      return
      endif
          
c--- CUT on the isolation of the leptons
      do i=1,2
      do j=7,maxparts ! jets correspond to indices 7,...,maxparts
        if (R(p,ilept(i),j) .lt. Risol) then
        if (pt(j,p) .gt. pt(ilept(i),p)*threshisol) then
          passed=.false.
          return
        endif
      endif
      enddo
      enddo
      
c--- print warning if there are not 2 neutrinos     
      if (inu .ne. 2) then
        write(6,*) 'WARNING: cdfhwwcuts.f must be upgraded, <2 neuts'
      stop
      endif     
      
c--- form the missing et vector 
      do j=1,2
        etvec(j)=etvec(j)+p(ineut(1),j)+p(ineut(2),j)
      enddo
      missinget=dsqrt(etvec(1)**2+etvec(2)**2)
      
c--- find the object (charged lepton or jet) that is closest in
c--- azimuthal angle to missing Et vector
      phimin=10d0
      do i=3,maxparts
        if ((i .ne. ineut(1)) .and. (i .ne. ineut(2))) then
          phi= (p(i,1)*etvec(1)+p(i,2)*etvec(2))
     .         /dsqrt((p(i,1)**2+p(i,2)**2)*missinget**2)
          if (phi .gt. +0.9999999D0) phi=+1D0
          if (phi .lt. -0.9999999D0) phi=-1D0
          phi=dacos(phi)
        endif
        if (phi .lt. phimin) phimin=phi
      enddo

c--- convert missing Et into missing Et "star" (or "spec")
      if (phi .lt. pi/2d0) then
        missinget=missinget*dsin(phi)
      endif

c--- CUT on missing Et 
      if (missinget .lt. missingetcut) then
        passed=.false.
      return
      endif

! OPTIONAL mt cut 
      m45=0d0
      do i=1,4
         if(i.ne.4) then 
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2 
         endif
      enddo


      mtmin=0d0*hmass
      mtmax=hmass 
      mt45=0d0 
      mt45=(dsqrt(dsqrt(pttwo(4,5,p)**2+m45)+etmiss(p,et_vec))**2)
!      write(6,*) mt45
 !     mt45=dsqrt(max(0d0,dsqrt(mt45-ptsq))) 
 !     write(6,*) mt45 
 !     pause 
!------ Mt cuts
      if((mt45.lt.mtmin).or.(mt45.gt.mtmax)) then 
         passed=.false. 
         return 
      endif



      return      
      end
      
     
