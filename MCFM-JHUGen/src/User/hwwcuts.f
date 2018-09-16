      subroutine hwwcuts(p,maxparts,dphi_ll,m_ll,mtrans,scut1,scut2)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic H->WW search cuts and then returns
c--- the values of dphi_ll,m_ll,mtrans,scut1,scut2 for the
c--- plotting routine; a point that fails the cuts returns mtrans=-1._dp

      include 'constants.f'
      include 'mxpart.f'
      include 'plabel.f'
      integer:: ilept(2),ic,j,maxparts
      real(dp):: pt,etarap,p(mxpart,4),etvec(2),plept(4),
     & missinget,dphi_ll,m_ll,mtrans,scut1,scut2,cosdphi,maxpt,minpt

c--- initialize all variables
      dphi_ll=-one
      m_ll=-one
      mtrans=-one
      scut1=-one
      scut2=-one
      etvec(1)=zip
      etvec(2)=zip

************************** START BASIC CUTS ****************************

      ic=0
c--- find two leptons that satisfy the pt and rapidity cuts
      do j=3,maxparts
        if (     (plabel(j)=='el') .or. (plabel(j)=='ea')
     &      .or. (plabel(j)=='ml') .or. (plabel(j)=='ma')) then
          if (       (pt(j,p) > 20._dp) .and.
     &      (abs(etarap(j,p)) < 2.5_dp)) then
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
        if ((j==1) .or. (j==2))
     &  etvec(j)=etvec(j)-p(ilept(1),j)-p(ilept(2),j)
      enddo

      missinget=sqrt(etvec(1)**2+etvec(2)**2)

c--- perform the cut on missing Et, rejecting less than 30 GeV
      if (missinget < 25._dp) goto 999

c---  perform veto if the event contains jets of 30 GeV with |eta|<3
      do j=3,maxparts
        if (     (plabel(j) == 'pp') .or. (plabel(j) == 'qj')
     &      .or. (plabel(j) == 'bq') .or. (plabel(j) == 'ba')) then
          if (       (pt(j,p) > 20._dp) .and.
     &       (abs(etarap(j,p)) < 3._dp)) then
            goto 999
          endif
        endif
      enddo

c--- opening angle between the leptons in the transverse plane
      dphi_ll=
     &   (p(ilept(1),1)*p(ilept(2),1)+p(ilept(1),2)*p(ilept(2),2))
     &   /sqrt((p(ilept(1),1)**2+p(ilept(1),2)**2)
     &         *(p(ilept(2),1)**2+p(ilept(2),2)**2))
      if (dphi_ll < -0.999999999_dp) dphi_ll=-one
      dphi_ll=acos(dphi_ll)

      if (dphi_ll > pi/four) then
        goto 999
      endif

c--- dilepton invariant mass
      m_ll=sqrt(plept(4)**2-plept(1)**2-plept(2)**2-plept(3)**2)

      if (m_ll > 35._dp) then
        goto 999
      endif

      maxpt=max(pt(ilept(1),p),pt(ilept(2),p))
      minpt=min(pt(ilept(1),p),pt(ilept(2),p))

      if (minpt < 25._dp) then
        goto 999
      endif

      if ((maxpt < 35._dp) .or. (maxpt > 50._dp)) then
        goto 999
      endif

*************************** END BASIC CUTS *****************************

******************** CALCULATE QUANTITIES TO PLOT **********************

c--- opening angle between the leptons in the transverse plane
      dphi_ll=
     &   (p(ilept(1),1)*p(ilept(2),1)+p(ilept(1),2)*p(ilept(2),2))
     &   /sqrt((p(ilept(1),1)**2+p(ilept(1),2)**2)
     &         *(p(ilept(2),1)**2+p(ilept(2),2)**2))
      if (dphi_ll < -0.999999999_dp) dphi_ll=-one
      dphi_ll=acos(dphi_ll)

c--- dilepton invariant mass
      m_ll=sqrt(plept(4)**2-plept(1)**2-plept(2)**2-plept(3)**2)

c--- transverse mass
c---    first compute the cosine of the azimuthal angle between the
c---    lepton system and the missing Et (often, this will simply be pi)
      cosdphi=
     &   (plept(1)*etvec(1)+plept(2)*etvec(2))
     &   /sqrt((plept(1)**2+plept(2)**2)
     &         *(etvec(1)**2+etvec(2)**2))
c---    transverse mass calculation
      mtrans=two*sqrt(plept(1)**2+plept(2)**2)*missinget*(one-cosdphi)
      mtrans=sqrt(max(mtrans,zip))

********************** END OF QUANTITIES TO PLOT ***********************

********************* START OF CUT CROSS SECTIONS **********************

c--- values of scut1,scut2 less than zero indicate that this point
c--- should not contribute to those cross sections

c--- require the azimuthal angle to be less than 45 degrees
      if (dphi_ll > pi/four) goto 999

c--- require the dilepton mass to be less than 35 GeV
      if (m_ll > 35._dp) goto 999

c--- scut1 corresponds to the cross section with 125 < mtrans < 155
      if ((mtrans >= 125._dp) .and. (mtrans <= 155._dp)) scut1=one

c--- scut2 corresponds to the cross section with 140 < mtrans < 180
      if ((mtrans >= 140._dp) .and. (mtrans <= 180._dp)) scut2=one

********************** END OF CUT CROSS SECTIONS ***********************

      return


c--- this is the alternate return if the cuts are failed
  999 continue

      return

      end


