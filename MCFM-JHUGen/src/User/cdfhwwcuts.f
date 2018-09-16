      subroutine cdfhwwcuts(p,maxparts,passed)
      implicit none
      include 'types.f'
c--- given the event momenta in p (with maxparts entries),
c--- performs H->WW search cuts along the lines of CDF Note 9887;
c--- a point that fails the cuts returns passed=.false.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'plabel.f'
      integer:: ilept(2),ile,ineut(2),inu,i,j,maxparts
      logical:: passed
      real(dp):: pt,etarap,p(mxpart,4),etvec(2),
     & missinget,m_ll,R,pt1,pt2,eta1,eta2,phi,phimin
c--- pt and eta thresholds for trigger lepton
      real(dp):: trigpt,trigeta
c--- pt and eta thresholds for second lepton (should be looser)
      real(dp):: scndpt,scndeta
c--- minimum threshold for dilepton invariant mass
      real(dp):: m_llcut
c--- minimum threshold for missing Et "star" or "spec" (see defn below)
      real(dp):: missingetcut
c--- size of cone for lepton isolation
      real(dp):: Risol
c--- fractional pt threshold for lepton isolation
      real(dp):: threshisol
      real(dp):: mtmax,mtmin,mt45,m45,et_vec(4),etmiss,pttwo
      data trigpt,trigeta/20._dp,0.8_dp/
      data scndpt,scndeta/10._dp,1.1_dp/
      data m_llcut/16._dp/
      data missingetcut/25._dp/
      data Risol/0.4_dp/
      data threshisol/0.1_dp/

c--- default behaviour is that the cuts have been passed
      passed=.true.

c--- initialize all variables
      etvec(1)=zip
      etvec(2)=zip

************************** START BASIC CUTS ****************************

      ile=0
      inu=0
c--- identify the two leptons and neutrinos
      do j=3,maxparts
        if (     (plabel(j)=='el') .or. (plabel(j)=='ea')
     &      .or. (plabel(j)=='ml') .or. (plabel(j)=='ma')) then
          ile=ile+1
          ilept(ile)=j
        elseif  ((plabel(j)=='nl') .or. (plabel(j)=='na')) then
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
      if ((  (pt1 > trigpt) .and. (eta1 < trigeta)
     &  .and.(pt2 > scndpt) .and. (eta2 < scndeta)) .or.
     &    (  (pt2 > trigpt) .and. (eta2 < trigeta)
     &  .and.(pt1 > scndpt) .and. (eta1 < scndeta)) ) then
        continue ! cut passed
      else
        passed=.false.
      return
      endif

c--- compute invariant mass of leptons
      m_ll=sqrt((p(ilept(1),4)+p(ilept(2),4))**2
     &          -(p(ilept(1),1)+p(ilept(2),1))**2
     &          -(p(ilept(1),2)+p(ilept(2),2))**2
     &          -(p(ilept(1),3)+p(ilept(2),3))**2)

c--- CUT on dilepton invariant mass
      if (m_ll < m_llcut) then
        passed=.false.
      return
      endif

c--- CUT on the isolation of the leptons
      do i=1,2
      do j=7,maxparts ! jets correspond to indices 7,...,maxparts
        if (R(p,ilept(i),j) < Risol) then
        if (pt(j,p) > pt(ilept(i),p)*threshisol) then
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
      missinget=sqrt(etvec(1)**2+etvec(2)**2)

c--- find the object (charged lepton or jet) that is closest in
c--- azimuthal angle to missing Et vector
      phimin=10._dp
      do i=3,maxparts
        if ((i .ne. ineut(1)) .and. (i .ne. ineut(2))) then
          phi= (p(i,1)*etvec(1)+p(i,2)*etvec(2))
     &         /sqrt((p(i,1)**2+p(i,2)**2)*missinget**2)
          if (phi > +0.9999999_dp) phi=+1._dp
          if (phi < -0.9999999_dp) phi=-1._dp
          phi=acos(phi)
        endif
        if (phi < phimin) phimin=phi
      enddo

c--- convert missing Et into missing Et "star" (or "spec")
      if (phi < pi/2._dp) then
        missinget=missinget*sin(phi)
      endif

c--- CUT on missing Et
      if (missinget < missingetcut) then
        passed=.false.
      return
      endif

! OPTIONAL mt cut
      m45=zip
      do i=1,4
         if(i.ne.4) then
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2
         endif
      enddo


      mtmin=zip*hmass
      mtmax=hmass
      mt45=zip
      mt45=(sqrt(sqrt(pttwo(4,5,p)**2+m45)+etmiss(p,et_vec))**2)
!      write(6,*) mt45
 !     mt45=sqrt(max(zip,sqrt(mt45-ptsq)))
 !     write(6,*) mt45
 !     pause
!------ Mt cuts
      if((mt45<mtmin).or.(mt45>mtmax)) then
         passed=.false.
         return
      endif



      return
      end


