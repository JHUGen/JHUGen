      subroutine hwwjetplots(eventpart,tag,p,wt,wt2)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with eventparts entries),
c--- makes some plots for the H->WW + jet search 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      include 'plabel.f'
      integer:: ilept(2),ic,ijet(3),jc,j,eventpart,n,nplotmax,i1,i2
      real(dp):: pt,etarap,p(mxpart,4),etvec(2),plept(4),
     & dphi_ll,m_ll,missinget,wt,wt2,etajet,etajet2,mtrans,
     & ptjet,ptjet2,Etll,Etmiss,ptH,pt1,pt2,pt3
      integer tag
      logical:: first
      common/nplotmax/nplotmax
      data first/.true./
      save first,ic,jc,ilept,ijet
ccccc!$omp threadprivate(/nplotmax/)
          
c--- initialize all variables  
      dphi_ll=-one
      m_ll=-one
      etajet=99._dp
      etvec(1)=zip
      etvec(2)=zip

      if (first) then
c--- if entering for the first time, identify positions of leptons and jets
        ic=0      
        jc=0
c--- find the two leptons and 1,2 or 3 jets 
        do j=3,eventpart
          if (     (plabel(j)=='el') .or. (plabel(j)=='ea')
     &        .or. (plabel(j)=='ml') .or. (plabel(j)=='ma')) then
              ic=ic+1
              if (ic <= 2) ilept(ic)=j
          endif
          if (     (plabel(j)=='bq') .or. (plabel(j)=='ba')
     &        .or. (plabel(j)=='pp')) then
              jc=jc+1
              if (jc <= 2) ijet(jc)=j
        endif
        enddo
c--- error if there are not 2 leptons      
        if (ic .ne. 2) then
          write(6,*) 'Expected to find 2 leptons to plot - found',ic
          stop
        endif     
c--- error if there are not 1 or 2 jets      
        if ((jc == 0) .or. (jc > 3)) then
          write(6,*) 'Expected to find 1,2 or 3 jets - found',jc
          stop
        endif     
      goto 99
      endif
      
******************** CALCULATE QUANTITIES TO PLOT **********************

c--- form the combined lepton momentum and the missing et vector 
      do j=1,4
        plept(j)=p(ilept(1),j)+p(ilept(2),j)
        if ((j==1) .or. (j==2)) then
        etvec(j)=etvec(j)-plept(j)-p(ijet(1),j)
        if (jc == 2) etvec(j)=etvec(j)-p(ijet(2),j)
      endif
      enddo

      missinget=sqrt(etvec(1)**2+etvec(2)**2)
      
c--- opening angle between the leptons in the transverse plane
      dphi_ll=
     &   (p(ilept(1),1)*p(ilept(2),1)+p(ilept(1),2)*p(ilept(2),2))
     &   /sqrt((p(ilept(1),1)**2+p(ilept(1),2)**2)
     &         *(p(ilept(2),1)**2+p(ilept(2),2)**2))
      if (dphi_ll < -0.999999999_dp) dphi_ll=-one
      dphi_ll=acos(dphi_ll) 
            
c--- dilepton invariant mass
      m_ll=sqrt(plept(4)**2-plept(1)**2-plept(2)**2-plept(3)**2)

c--- rapidities and pts of leading (highest-pt) jet and sub-leading one
      etajet2=99._dp
      ptjet2=-one
      if     (jets == 1) then
        etajet=etarap(ijet(1),p)
        ptjet=pt(ijet(1),p)
      elseif (jets == 2) then
        ptjet=pt(ijet(1),p)
      ptjet2=pt(ijet(2),p)
        if (ptjet > ptjet2) then
          etajet=etarap(ijet(1),p)
          etajet2=etarap(ijet(2),p)
        else
        etajet=ptjet
        ptjet=ptjet2
        ptjet2=etajet
        etajet=etarap(ijet(2),p)
        etajet2=etarap(ijet(1),p)
        endif
      elseif (jets == 3) then
c--- sort for 3 jets 
        pt1=pt(ijet(1),p)
        pt2=pt(ijet(2),p)
        pt3=pt(ijet(3),p)
        if ((pt1 > pt2) .and. (pt1 > pt3)) then
          i1=1
          if (pt2 > pt3) then
            i2=2
c            i3=3
          else
            i2=3
c            i3=2
          endif
        endif
        if ((pt2 > pt1) .and. (pt2 > pt3)) then
          i1=2
          if (pt1 > pt3) then
            i2=1
c            i3=3
          else
            i2=3
c            i3=1
          endif
        endif
        if ((pt3 > pt1) .and. (pt3 > pt2)) then
           i1=3
          if (pt1 > pt2) then
            i2=1
c            i3=2
          else
            i2=2
c            i3=1
          endif
        endif
        ptjet=pt(ijet(i1),p)
        ptjet2=pt(ijet(i2),p)
        etajet=etarap(ijet(i1),p)
        etajet2=etarap(ijet(i2),p)
      endif
      
c--- NOT USED ANY MORE (SEE BELOW)
cc--- transverse mass
cc---    first compute the cosine of the azimuthal angle between the
cc---    lepton system and the missing Et (often, this will simply be pi)
c      cosdphi=
c     &   (plept(1)*etvec(1)+plept(2)*etvec(2))
c     &   /sqrt((plept(1)**2+plept(2)**2)
c     &         *(etvec(1)**2+etvec(2)**2))
cc---    transverse mass calculation
c      mtrans=2._dp*sqrt(plept(1)**2+plept(2)**2)*missinget*(one-cosdphi)
c      mtrans=sqrt(max(mtrans,zip))

c--- Transverse mass of Klamke and Zeppenfeld,
c---  c.f. Eqs. 13 and 14 of hep-ph/0703202.
      Etll=sqrt(
     & +(p(ilept(1),4)+p(ilept(2),4))**2
     & -(p(ilept(1),3)+p(ilept(2),3))**2)
      Etmiss=sqrt(
     & +etvec(1)**2+etvec(2)**2
     & +(p(ilept(1),4)+p(ilept(2),4))**2
     & -(p(ilept(1),1)+p(ilept(2),1))**2
     & -(p(ilept(1),2)+p(ilept(2),2))**2
     & -(p(ilept(1),3)+p(ilept(2),3))**2)
      mtrans=sqrt((Etmiss+Etll)**2
     & -(p(ilept(1),1)+p(ilept(2),1)+etvec(1))**2
     & -(p(ilept(1),2)+p(ilept(2),2)+etvec(2))**2)

      ptH=sqrt((plept(1)+etvec(1))**2+(plept(2)+etvec(2))**2)

********************** END OF QUANTITIES TO PLOT ***********************

   99 continue

      n=1
      call bookplot(n,tag,'eta jet1',etajet,wt,wt2,
     &              -4.5_dp,4.5_dp,0.3_dp,'lin')
      n=n+1
      call bookplot(n,tag,'dphi_ll',dphi_ll,wt,wt2,
     &              zip,3.1416_dp,0.10472_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m_ll',m_ll,wt,wt2,
     &              zip,300._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'mtrans',mtrans,wt,wt2,
     &              zip,1000._dp,33.333_dp,'log')
      n=n+1
      call bookplot(n,tag,'Et miss',missinget,wt,wt2,
     &              zip,1000._dp,33.333_dp,'log')
      n=n+1
      call bookplot(n,tag,'pt jet1',ptjet,wt,wt2,
     &              30._dp,1000._dp,33.333_dp,'log')
      n=n+1
      call bookplot(n,tag,'pt Higgs',ptH,wt,wt2,
     &              zip,1000._dp,33.333_dp,'log')
      n=n+1
      call bookplot(n,tag,'eta jet2',etajet2,wt,wt2,
     &              -4.5_dp,4.5_dp,0.3_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt jet2',ptjet2,wt,wt2,
     &              30._dp,1000._dp,33.333_dp,'log')
      n=n+1
      
      n=n-1
      
c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return      
      end
      
     
