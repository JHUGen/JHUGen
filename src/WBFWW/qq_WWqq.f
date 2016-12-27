      subroutine qq_WWqq(p,msq)
      implicit none
c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W(p3,p4)+W(p5,p6)+q(p7)+q(p8);
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      !include 'first.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'pid_pdg.f'
      include 'WWbits.f'
      integer nmax,jmax
      parameter(jmax=12,nmax=10)
      integer j,k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq,
     & i3,i4,i5,i6,nfinc
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,mult,
     & colfac34_56
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2)
      logical doHO,doBO,comb1278ok
      logical isALepton,isANeutrino
      parameter(spinavge=0.25d0,stat=0.5d0,nfinc=4)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO,mult

      msq(:,:)=0d0

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      cwmass2=dcmplx(wmass**2,0d0)
      czmass2=dcmplx(zmass**2,0d0)
      cxw=dcmplx(xw,0d0)

      doHO=.false.
      doBO=.false.
      if     (runstring(4:5) .eq. 'HO') then
        doHO=.true.
      elseif (runstring(4:5) .eq. 'BO') then
        doBO=.true.
      endif
      if (doHO) then
        Hbit=cone
        Bbit=czip
      elseif (doBO) then
        Hbit=czip
        Bbit=cone
      else
        Hbit=cone
        Bbit=cone
      endif

c--- rescaling factor for Higgs amplitudes, if anomalous Higgs width
      mult=1d0
      if (anom_Higgs) then
         mult=chi_higgs**2
      endif
      Hbit=mult*Hbit

C---call plabel/pdgid conversion
      call convertPLabelsToPDGIds()

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

c---color factors for W decays
      colfac34_56=1d0
      if (abs(pid_pdg(3)).ge.0 .and. abs(pid_pdg(3)).le.5) then
        colfac34_56=colfac34_56*xn
      endif
      if (abs(pid_pdg(5)).ge.0 .and. abs(pid_pdg(5)).le.5) then
        colfac34_56=colfac34_56*xn
      endif

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:)=czip
      ampa(:,:,:)=czip
      ampb(:,:,:)=czip

C--   MARKUS: adding switches to remove VH or VBF contributions
      ! No VH-like diagram
      !if( (vvhvvtoggle_vbfvh.eq.0) .and. (j.ge.9) ) cycle
      ! No VBF-like diagram
      !if( (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.8) ) cycle
      if( (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.4) ) cycle
      ! U. Sarica: Test the combination
      call testWBFVVApartComb(j1(j),j2(j),j7(j),j8(j),comb1278ok)
      if (.not.comb1278ok) cycle
      if( (
     &    (isALepton(abs(pid_pdg(7))) .or. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (isALepton(abs(pid_pdg(8))) .or. isANeutrino(abs(pid_pdg(8))))
     &    ) .and. j.lt.9
     & ) then
         cycle
      endif


c---  Call the VVWW amplitudes
c---  Note 3654, this is the special ordering to agree with Madgraph
      ! This orders the decay bosons as W-W+ as in process.DAT
      call getVVWWamps(amp,ampa,ampb,p,za,zb,zab,zba,
     & j1(j),j2(j),3,6,5,4,j7(j),j8(j),doHO,doBO)

      ! Kill amp or ampa in j=5,6,7,8 for VH
      ! Kill amp or ampa in j=9,10,11,12 for VBF
      if( (
     & (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.8)
     & ) .or. (
     & (vvhvvtoggle_vbfvh.eq.0) .and. (j.ge.9)
     & )
     &  ) then
         amp(:,:,:)=czip
         ampa(:,:,:)=czip
      ! Kill ampb in j=9,10,11,12 for VH
      ! Kill ampb in j=5,6,7,8 for VBF
      else if( (
     & (vvhvvtoggle_vbfvh.eq.1) .and. (j.ge.9)
     & ) .or. (
     & (vvhvvtoggle_vbfvh.eq.0) .and. (j.le.8 .and. j.ge.5)
     & )
     &  ) then
         ampb(:,:,:)=czip
      endif
      if(
     &    (isALepton(abs(pid_pdg(7))) .or. isANeutrino(abs(pid_pdg(7))))
     &    .or.
     &    (isALepton(abs(pid_pdg(8))) .or. isANeutrino(abs(pid_pdg(8))))
     & ) then
         ampa(:,:,:)=czip
         ampb(:,:,:)=czip
      endif

C-----setup for (dqcq_dqcq)
      do h1=1,2
      do h2=1,2
      temp(1,4)=temp(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_dqcq,h1,h2)
     & *dconjg(amp(dqcq_dqcq,h1,h2)))
      enddo
      enddo

C-----setup for (dqcq_uqsq)
      tempw(1,4)=tempw(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_uqsq,1,1)
     & *dconjg(amp(dqcq_uqsq,1,1)))

C-----setup for (uqcq_uqcq)
      do h1=1,2
      do h2=1,2
      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *dble(amp(uqcq_uqcq,h1,h2)
     & *dconjg(amp(uqcq_uqcq,h1,h2)))
      enddo
      enddo


C-----setup for (dqsq_dqsq)
      do h1=1,2
      do h2=1,2
      temp(1,3)=temp(1,3)+esq**6*spinavge
     &   *dble(amp(dqsq_dqsq,h1,h2)
     & *dconjg(amp(dqsq_dqsq,h1,h2)))
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

C-----setup for (dqdq_dqdq)
      do h1=1,2
      do h2=1,2
      temp(1,1)=temp(1,1)+esq**6*spinavge
     &   *dble(ampa(dqdq_dqdq,h1,h2)
     & *dconjg(ampa(dqdq_dqdq,h1,h2)))
      temp(1,1)=temp(1,1)+esq**6*spinavge
     &   *dble(ampb(dqdq_dqdq,h1,h2)
     & *dconjg(ampb(dqdq_dqdq,h1,h2)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*esq**6*spinavge
     &   *dble(ampa(dqdq_dqdq,h1,h2)
     & *dconjg(ampb(dqdq_dqdq,h1,h2)))
      endif
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)

C-----setup for (uquq_uquq)
      do h1=1,2
      do h2=1,2
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *dble(ampa(uquq_uquq,h1,h2)
     & *dconjg(ampa(uquq_uquq,h1,h2)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *dble(ampb(uquq_uquq,h1,h2)
     & *dconjg(ampb(uquq_uquq,h1,h2)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*esq**6*spinavge
     &   *dble(ampa(uquq_uquq,h1,h2)
     & *dconjg(ampb(uquq_uquq,h1,h2)))
      endif
      enddo
      enddo
      temp(4,4)=temp(2,2)

      do h1=1,2
      do h2=1,2
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2)
     & *dconjg(ampa(dquq_dquq,h1,h2)))
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampb(dquq_dquq,h1,h2)
     & *dconjg(ampb(dquq_dquq,h1,h2)))
      if ((h1 .eq. 1) .and. (h2 .eq. 1)) then
      temp(1,2)=temp(1,2)-2d0/xn*esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2)
     & *dconjg(ampb(dquq_dquq,h1,h2)))
      endif
      enddo
      enddo
      temp(3,4)=temp(1,2)

C-----setup for (uqbq_uqbq)
      do h1=1,2
      do h2=1,2
      temp(2,5)=temp(2,5)+esq**6*spinavge
     &   *dble(amp(uqbq_uqbq,h1,h2)
     & *dconjg(amp(uqbq_uqbq,h1,h2)))
      enddo
      enddo
      temp(4,5)=temp(2,5)
      temp(2,3)=temp(2,5)

C-----setup for (uqsq_dqcq)
      tempw(2,3)=tempw(2,3)+esq**6*spinavge
     &   *dble(amp(uqsq_dqcq,1,1)
     & *dconjg(amp(uqsq_dqcq,1,1)))

      do k=1,nf;do l=1,nf
      if (
     & .not.(
     &    (
     &    (pid_pdg(j1(j)).eq.0)
     &    .or. (
     &    j1(j).le.2 .and. (
     &    pid_pdg(j1(j)).eq.k
     &    .or.
     &    (
     &    modulo(k,2).eq.0 .and.
     &    (pid_pdg(j1(j)).eq.12
     &    .or. pid_pdg(j1(j)).eq.14
     &    .or. pid_pdg(j1(j)).eq.16)
     &    )
     &    .or.
     &    (
     &    modulo(k,2).eq.1 .and.
     &    (pid_pdg(j1(j)).eq.11
     &    .or. pid_pdg(j1(j)).eq.13
     &    .or. pid_pdg(j1(j)).eq.15)
     &    )
     &    )
     &    )
     &    .or. (
     &    j1(j).ge.7 .and. (
     &    pid_pdg(j1(j)).eq.-k
     &    .or.
     &    (
     &    modulo(k,2).eq.0 .and.
     &    (pid_pdg(j1(j)).eq.-12
     &    .or. pid_pdg(j1(j)).eq.-14
     &    .or. pid_pdg(j1(j)).eq.-16)
     &    )
     &    .or.
     &    (
     &    modulo(k,2).eq.1 .and.
     &    (pid_pdg(j1(j)).eq.-11
     &    .or. pid_pdg(j1(j)).eq.-13
     &    .or. pid_pdg(j1(j)).eq.-15)
     &    )
     &    )
     &    )
     &    ) .and. (
     &    (pid_pdg(j2(j)).eq.0)
     &    .or. (
     &    j2(j).le.2 .and. (
     &    pid_pdg(j2(j)).eq.l
     &    .or.
     &    (
     &    modulo(l,2).eq.0 .and.
     &    (pid_pdg(j2(j)).eq.12
     &    .or. pid_pdg(j2(j)).eq.14
     &    .or. pid_pdg(j2(j)).eq.16)
     &    )
     &    .or.
     &    (
     &    modulo(l,2).eq.1 .and.
     &    (pid_pdg(j2(j)).eq.11
     &    .or. pid_pdg(j2(j)).eq.13
     &    .or. pid_pdg(j2(j)).eq.15)
     &    )
     &    )
     &    )
     &    .or. (
     &    j2(j).ge.7 .and. (
     &    pid_pdg(j2(j)).eq.-l
     &    .or.
     &    (
     &    modulo(l,2).eq.0 .and.
     &    (pid_pdg(j2(j)).eq.-12
     &    .or. pid_pdg(j2(j)).eq.-14
     &    .or. pid_pdg(j2(j)).eq.-16)
     &    )
     &    .or.
     &    (
     &    modulo(l,2).eq.1 .and.
     &    (pid_pdg(j2(j)).eq.-11
     &    .or. pid_pdg(j2(j)).eq.-13
     &    .or. pid_pdg(j2(j)).eq.-15)
     &    )
     &    )
     &    )
     &    )
     & )
     & ) then
         temp(k,l)=zip
         tempw(k,l)=zip
      endif
      enddo;enddo
      if (
     & isANeutrino(abs(pid_pdg(7))) .and.
     & abs(pid_pdg(8)).eq.abs(pid_pdg(7))
     & ) then
         tempw(:,:)=zip
         if(j.eq.10 .or. j.eq.12) then
            temp(4,5)=zip
         endif
      else if (
     & isALepton(abs(pid_pdg(7))) .and.
     & abs(pid_pdg(8)).eq.abs(pid_pdg(7))
     & ) then
         tempw(:,:)=zip
         if (j.eq.9 .or. j.eq.11) then
            temp(1:2,3)=zip
         else if (j.eq.10 .or. j.eq.12) then
            temp(1,3:5)=zip
         endif
      else if (
     & (isALepton(abs(pid_pdg(7))) .and.
     & isANeutrino(abs(pid_pdg(8)))) .or.
     & (isALepton(abs(pid_pdg(8))) .and.
     & isANeutrino(abs(pid_pdg(7))))
     & ) then
         temp(:,:)=zip
      endif

c---- Multiply by the decay color factor
      temp(:,:) = temp(:,:)*colfac34_56
      tempw(:,:) = tempw(:,:)*colfac34_56

c--- fill matrix elements
      if (j.eq.1) then
      do k=1,nf
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nf
      msq(k,l)=temp(k,l)
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         msq(2,3)=msq(2,3)+tempw(2,3)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         msq(1,4)=msq(1,4)+tempw(1,4)
      endif

      elseif (j.eq.2) then
      do k=1,nf
      do l=k+1,nf
      msq(l,k)=temp(k,l)
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         msq(3,2)=msq(3,2)+tempw(2,3)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         msq(4,1)=msq(4,1)+tempw(1,4)
      endif

      elseif (j.eq.3) then
      do k=-nf,-1
      msq(k,k)=temp(-k,-k)*stat
      do l=k+1,-1
      msq(k,l)=temp(-l,-k)
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         msq(-3,-2)=msq(-3,-2)+tempw(1,4)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         msq(-4,-1)=msq(-4,-1)+tempw(2,3)
      endif

      elseif (j.eq.4) then
      do k=-nf,-1
      do l=k+1,-1
      msq(l,k)=temp(-l,-k)
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     v .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         msq(-2,-3)=msq(-2,-3)+tempw(1,4)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         msq(-1,-4)=msq(-1,-4)+tempw(2,3)
      endif

c--- qbar-q
      elseif (j.eq.5) then
      do k=-nf,-1
      msq(k,-k)=temp(-k,-k)
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      msq(k,l)=temp(-k,l)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         msq(-1,3)=msq(-1,3)+tempw(2,3)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         msq(-2,4)=msq(-2,4)+tempw(1,4)
      endif

c--- qbar-q
      elseif (j.eq.6) then
      do k=-nf,-1
      do l=1,nf
      if (abs(k) .gt. abs(l)) then
      msq(k,l)=temp(l,-k)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         msq(-3,1)=msq(-3,1)+tempw(1,4)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         msq(-4,2)=msq(-4,2)+tempw(2,3)
      endif

c--- q-qbar
      elseif (j.eq.7) then
      do k=-nf,-1
      msq(-k,k)=temp(-k,-k)
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      msq(l,k)=temp(-k,l)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         msq(3,-1)=msq(3,-1)+tempw(2,3)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & ) then
         msq(4,-2)=msq(4,-2)+tempw(1,4)
      endif

c--- q-qbar
      elseif (j.eq.8) then
      do k=-nf,-1
      do l=-nf,-1
      if (abs(k) .lt. abs(l)) then
      msq(-k,l)=temp(-k,-l)
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         msq(1,-3)=msq(1,-3)+tempw(1,4)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & ) then
         msq(2,-4)=msq(2,-4)+tempw(2,3)
      endif

c--- q-qbar extra pieces
      elseif (j.eq.9) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      msq(k,-k)=msq(k,-k)+temp(k,l)
      endif
      enddo
      enddo
      msq(3,-4)=msq(1,-2)
      msq(4,-3)=msq(2,-1)
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         msq(1,-2)=msq(1,-2)+tempw(1,4) ! d u~ -> c~ s
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         msq(3,-4)=msq(3,-4)+tempw(1,4) ! s c~ -> u~ d
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         msq(2,-1)=msq(2,-1)+tempw(2,3) ! u d~ -> s~ c
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         msq(4,-3)=msq(4,-3)+tempw(2,3) ! c s~ -> d~ u
      endif

c--- q-qbar extra pieces
      elseif (j.eq.10) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      msq(k,-k)=msq(k,-k)+temp(l,k)
      endif
      enddo
      enddo

c--- qbar-q extra pieces
      elseif (j.eq.11) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      msq(-k,k)=msq(-k,k)+temp(k,l)
      endif
      enddo
      enddo
      msq(-4,3)=msq(-2,1)
      msq(-3,4)=msq(-1,2)
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         msq(-2,1)=msq(-2,1)+tempw(1,4) ! u~ d -> c~ s
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2
     &                  .or. isANeutrino(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1
     &                  .or. isALepton(pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1
     &                  .or. isALepton(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2
     &                  .or. isANeutrino(-pid_pdg(8)))
     & )
     & ) then
         msq(-4,3)=msq(-4,3)+tempw(1,4) ! c~ s -> u~ d
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         msq(-1,2)=msq(-1,2)+tempw(2,3) ! d~ u -> s~ c
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2
     &                  .or. isANeutrino(pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1
     &                  .or. isALepton(-pid_pdg(8)))
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1
     &                  .or. isALepton(-pid_pdg(7))) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2
     &                  .or. isANeutrino(pid_pdg(8)))
     & )
     & ) then
         msq(-3,4)=msq(-3,4)+tempw(2,3) ! s~ c -> d~ u
      endif

c--- qbar-q extra pieces
      elseif (j.eq.12) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      msq(-k,k)=msq(-k,k)+temp(l,k)
      endif
      enddo
      enddo

      endif

      enddo

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end

