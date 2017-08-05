      subroutine qq_WWqqstrong(p,msq)
      implicit none
c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W(p3,p4)+W(p5,p6)+q(p7)+q(p8);
c--- with the t-channel exchange of a gluon.
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'pid_pdg.f'
      integer nmax,jmax
      parameter(jmax=12,nmax=10)
      integer j,k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq,
     & i3,i4,i5,i6
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,Colorfac,colfac34_56,ampsqfac,
     & msqgg(2)
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2)
      logical comb1278ok
      parameter(spinavge=0.25d0,stat=0.5d0,Colorfac=V/4d0/xn**2)
c--- note that this is the special ordering to agree with Madgraph
      parameter(i3=3,i4=6,i5=5,i6=4)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)

      msq(:,:)=0d0
      ampsqfac = esq**4*gsq**2*Colorfac*spinavge

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      cwmass2=dcmplx(wmass**2,0d0)
      czmass2=dcmplx(zmass**2,0d0)
      cxw=dcmplx(xw,0d0)

C---call plabel/pdgid conversion
      call convertPLabelsToPDGIds()

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

c---setup Z/A couplings from PDG ids
      call couplzajk()

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

      ! Test combination for qq'->q''q''' states
      call testWBFVVApartComb(j1(j),j2(j),j7(j),j8(j),comb1278ok)
      if (.not.comb1278ok) cycle

c--   Call the VVZZ amplitudes
      call getQQWWQQstrongamps(amp,ampa,ampb,za,zb,zab,zba,
     & j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j))

C-----setup for (dqcq_dqcq)
      do h1=1,2
      do h2=1,2
      temp(1,4)=temp(1,4)+ampsqfac
     &   *dble(amp(dqcq_dqcq,h1,h2)
     & *dconjg(amp(dqcq_dqcq,h1,h2)))
      enddo
      enddo

C-----setup for (dqcq_uqsq)
      tempw(1,4)=tempw(1,4)+ampsqfac
     &   *dble(amp(dqcq_uqsq,1,1)
     & *dconjg(amp(dqcq_uqsq,1,1)))

C-----setup for (uqcq_uqcq)
      do h1=1,2
      do h2=1,2
      temp(2,4)=temp(2,4)+ampsqfac
     &   *dble(amp(uqcq_uqcq,h1,h2)
     & *dconjg(amp(uqcq_uqcq,h1,h2)))
      enddo
      enddo


C-----setup for (dqsq_dqsq)
      do h1=1,2
      do h2=1,2
      temp(1,3)=temp(1,3)+ampsqfac
     &   *dble(amp(dqsq_dqsq,h1,h2)
     & *dconjg(amp(dqsq_dqsq,h1,h2)))
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

C-----setup for (dqdq_dqdq)
      do h1=1,2
      do h2=1,2
      temp(1,1)=temp(1,1)+ampsqfac
     &   *dble(ampa(dqdq_dqdq,h1,h2)
     & *dconjg(ampa(dqdq_dqdq,h1,h2)))
      temp(1,1)=temp(1,1)+ampsqfac
     &   *dble(ampb(dqdq_dqdq,h1,h2)
     & *dconjg(ampb(dqdq_dqdq,h1,h2)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*ampsqfac
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
      temp(2,2)=temp(2,2)+ampsqfac
     &   *dble(ampa(uquq_uquq,h1,h2)
     & *dconjg(ampa(uquq_uquq,h1,h2)))
      temp(2,2)=temp(2,2)+ampsqfac
     &   *dble(ampb(uquq_uquq,h1,h2)
     & *dconjg(ampb(uquq_uquq,h1,h2)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*ampsqfac
     &   *dble(ampa(uquq_uquq,h1,h2)
     & *dconjg(ampb(uquq_uquq,h1,h2)))
      endif
      enddo
      enddo
      temp(4,4)=temp(2,2)

C-----setup for (dquq_dquq)
      do h1=1,2
      do h2=1,2
      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampa(dquq_dquq,h1,h2)
     & *dconjg(ampa(dquq_dquq,h1,h2)))
      temp(1,2)=temp(1,2)+ampsqfac
     &   *dble(ampb(dquq_dquq,h1,h2)
     & *dconjg(ampb(dquq_dquq,h1,h2)))
      if (h1 .eq. h2) then
      temp(1,2)=temp(1,2)-2d0/xn*ampsqfac
     &   *dble(ampa(dquq_dquq,h1,h2)
     & *dconjg(ampb(dquq_dquq,h1,h2)))
      endif
      enddo
      enddo
      temp(3,4)=temp(1,2)

C-----setup for (uqbq_uqbq)
      do h1=1,2
      do h2=1,2
      temp(2,5)=temp(2,5)+ampsqfac
     &   *dble(amp(uqbq_uqbq,h1,h2)
     & *dconjg(amp(uqbq_uqbq,h1,h2)))
      enddo
      enddo
      temp(4,5)=temp(2,5)
      temp(2,3)=temp(2,5)

C-----setup for (uqsq_dqcq)
      tempw(2,3)=tempw(2,3)+ampsqfac
     &   *dble(amp(uqsq_dqcq,1,1)
     & *dconjg(amp(uqsq_dqcq,1,1)))

      do k=1,nf;do l=1,nf
      if (.not.(
     & (
     & (pid_pdg(j1(j)).eq.0)
     & .or. (j1(j).le.2 .and. pid_pdg(j1(j)).eq.k)
     & .or. (j1(j).ge.7 .and. pid_pdg(j1(j)).eq.-k)
     & ) .and. (
     & (pid_pdg(j2(j)).eq.0)
     & .or. (j2(j).le.2 .and. pid_pdg(j2(j)).eq.l)
     & .or. (j2(j).ge.7 .and. pid_pdg(j2(j)).eq.-l)
     & )
     & )) then
         temp(k,l)=zip
      endif
      enddo;enddo

c---- Multiply by the decay color factor
      temp(:,:) = temp(:,:)*colfac34_56
      tempw(:,:) = tempw(:,:)*colfac34_56

c--- fill matrix elements
      if (j.eq.1) then
      do k=1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,k,k,k,j,stat)
      endif
      do l=k+1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,k,l,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,2,3,2,3,j,1d0)
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
         call addtempwtomsq(msq,tempw,1,4,1,4,j,1d0)
      endif

      elseif (j.eq.2) then
      do k=1,nf
      do l=k+1,nf
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,k,l,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,3,2,2,3,j,1d0)
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
         call addtempwtomsq(msq,tempw,4,1,1,4,j,1d0)
      endif

      elseif (j.eq.3) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,k,-k,-k,j,stat)
      endif
      do l=k+1,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,-l,-k,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,-3,-2,1,4,j,1d0)
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
         call addtempwtomsq(msq,tempw,-4,-1,2,3,j,1d0)
      endif

      elseif (j.eq.4) then
      do k=-nf,-1
      do l=k+1,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,-l,-k,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,-2,-3,1,4,j,1d0)
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
         call addtempwtomsq(msq,tempw,-1,-4,2,3,j,1d0)
      endif

c--- qbar-q
      elseif (j.eq.5) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & ) .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,-k,-k,j,1d0)
      endif
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,-k,l,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,-1,3,2,3,j,1d0)
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
         call addtempwtomsq(msq,tempw,-2,4,1,4,j,1d0)
      endif

c--- qbar-q
      elseif (j.eq.6) then
      do k=-nf,-1
      do l=1,nf
      if (abs(k) .gt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,l,l,-k,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,-3,1,1,4,j,1d0)
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
         call addtempwtomsq(msq,tempw,-4,2,2,3,j,1d0)
      endif

c--- q-qbar
      elseif (j.eq.7) then
      do k=-nf,-1
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & ) .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.k)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,-k,-k,j,1d0)
      endif
      do l=1,nf
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.k) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,l,k,-k,l,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,3,-1,2,3,j,1d0)
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
         call addtempwtomsq(msq,tempw,4,-2,1,4,j,1d0)
      endif

c--- q-qbar
      elseif (j.eq.8) then
      do k=-nf,-1
      do l=-nf,-1
      if (abs(k) .lt. abs(l)) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-k)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,l,-k,-l,j,1d0)
      endif
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
         call addtempwtomsq(msq,tempw,1,-3,1,4,j,1d0)
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
         call addtempwtomsq(msq,tempw,2,-4,2,3,j,1d0)
      endif

c--- q-qbar extra pieces
      elseif (j.eq.9) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,1,-2,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,3,-4,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,2,-1,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,4,-3,2,3,j,1d0)
      endif

c--- q-qbar extra pieces
      elseif (j.eq.10) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,k,-k,l,k,j,1d0)
      endif
      endif
      enddo
      enddo

c--- qbar-q extra pieces
      elseif (j.eq.11) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,k,l,j,1d0)
      endif
      endif
      enddo
      enddo
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.3)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-2,1,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-2)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-4,3,1,4,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.4) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-3)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-3) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.4)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-1,2,2,3,j,1d0)
      endif
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.2) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-1)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-1) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.2)
     & )
     & ) then
         call addtempwtomsq(msq,tempw,-3,4,2,3,j,1d0)
      endif

c--- qbar-q extra pieces
      elseif (j.eq.12) then
      do k=1,nf
      do l=1,nf
      if (k .gt. l) then
      if (
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.-l)
     & )
     & .or.
     & (
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.-l) .and.
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.l)
     & )
     & ) then
         call addtemptomsq(msq,temp,-k,k,l,k,j,1d0)
      endif
      endif
      enddo
      enddo

      endif

      enddo

c--- 2-gluon amplitudes
      if (
     & ((pid_pdg(7).eq.0) .or. (pid_pdg(7).eq.21)) .and.
     & ((pid_pdg(8).eq.0) .or. (pid_pdg(8).eq.21))
     & ) then ! qqb/qbq->gg
         call qq2l2nuggamp(1,2,i3,i4,i5,i6,7,8,za,zb,msqgg)
         msq(1,-1)=msq(1,-1)+stat*aveqq*msqgg(1)
         msq(2,-2)=msq(2,-2)+stat*aveqq*msqgg(2)
         msq(3,-3)=msq(3,-3)+stat*aveqq*msqgg(1)
         msq(4,-4)=msq(4,-4)+stat*aveqq*msqgg(2)
         msq(5,-5)=msq(5,-5)+stat*aveqq*msqgg(1)

         call qq2l2nuggamp(2,1,i3,i4,i5,i6,7,8,za,zb,msqgg)
         msq(-1,1)=msq(-1,1)+stat*aveqq*msqgg(1)
         msq(-2,2)=msq(-2,2)+stat*aveqq*msqgg(2)
         msq(-3,3)=msq(-3,3)+stat*aveqq*msqgg(1)
         msq(-4,4)=msq(-4,4)+stat*aveqq*msqgg(2)
         msq(-5,5)=msq(-5,5)+stat*aveqq*msqgg(1)
      endif

      if (
     & ((pid_pdg(7).eq.0).and.(pid_pdg(8).le.5 .and. pid_pdg(8).ge.0))
     & .or.
     & ((pid_pdg(8).eq.0).and.(pid_pdg(7).ge.-5 .and. pid_pdg(7).le.0))
     & .or.
     & ((pid_pdg(8).le.5 .and. pid_pdg(8).ge.0)
     & .and. (pid_pdg(7).eq.(-pid_pdg(8))))
     & ) then ! gg->qbq
         call qq2l2nuggamp(7,8,i3,i4,i5,i6,1,2,za,zb,msqgg)
         do k=1,5,2 ! db-d,sb-s,bb-b
            if(
     &       ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) .and.
     &       ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k))
     &       ) then
               msq(0,0)=msq(0,0)+avegg*msqgg(1)
            endif
         enddo
         do k=2,4,2 ! ub-u,cb-c
            if(
     &       ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) .and.
     &       ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k))
     &       ) then
               msq(0,0)=msq(0,0)+avegg*msqgg(2)
            endif
         enddo
      endif

      if(
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.21) .and.
     & (pid_pdg(7).ge.-5 .and. pid_pdg(7).le.0)
     & ) then ! gqb->qbg / qbg->qbg
         ! gqb->qbg
         call qq2l2nuggamp(7,2,i3,i4,i5,i6,1,8,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(0,-k)=msq(0,-k)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(0,-k)=msq(0,-k)+aveqg*msqgg(2)
         endif
         enddo

         ! qbg->qbg
         call qq2l2nuggamp(7,1,i3,i4,i5,i6,2,8,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(-k,0)=msq(-k,0)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(-k,0)=msq(-k,0)+aveqg*msqgg(2)
         endif
         enddo
      else if(
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.21) .and.
     & (pid_pdg(8).ge.-5 .and. pid_pdg(8).le.0)
     & ) then ! gqb->gqb / qbg->gqb
         ! gqb->gqb
         call qq2l2nuggamp(8,2,i3,i4,i5,i6,1,7,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(0,-k)=msq(0,-k)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(0,-k)=msq(0,-k)+aveqg*msqgg(2)
         endif
         enddo

         ! qbg->gqb
         call qq2l2nuggamp(8,1,i3,i4,i5,i6,2,7,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(-k,0)=msq(-k,0)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(-k,0)=msq(-k,0)+aveqg*msqgg(2)
         endif
         enddo
      endif

      if(
     & (pid_pdg(7).eq.0 .or. pid_pdg(7).eq.21) .and.
     & (pid_pdg(8).le.5 .and. pid_pdg(8).ge.0)
     & ) then ! gq->gq / qg->gq
         ! gq->gq
         call qq2l2nuggamp(2,8,i3,i4,i5,i6,1,7,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(0,k)=msq(0,k)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(0,k)=msq(0,k)+aveqg*msqgg(2)
         endif
         enddo

         ! qg->gq
         call qq2l2nuggamp(1,8,i3,i4,i5,i6,2,7,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(k,0)=msq(k,0)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(8).eq.0) .or. (abs(pid_pdg(8)).eq.k)) then
            msq(k,0)=msq(k,0)+aveqg*msqgg(2)
         endif
         enddo
      else if(
     & (pid_pdg(8).eq.0 .or. pid_pdg(8).eq.21) .and.
     & (pid_pdg(7).le.5 .and. pid_pdg(7).ge.0)
     & ) then ! gq->qg / qg->qg
         ! gq->qg
         call qq2l2nuggamp(2,7,i3,i4,i5,i6,1,8,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(0,k)=msq(0,k)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(0,k)=msq(0,k)+aveqg*msqgg(2)
         endif
         enddo

         ! qg->qg
         call qq2l2nuggamp(1,7,i3,i4,i5,i6,2,8,za,zb,msqgg)
         do k=1,5,2 ! q=d,s,b
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(k,0)=msq(k,0)+aveqg*msqgg(1)
         endif
         enddo
         do k=2,4,2 ! q=u,c
         if ((pid_pdg(7).eq.0) .or. (abs(pid_pdg(7)).eq.k)) then
            msq(k,0)=msq(k,0)+aveqg*msqgg(2)
         endif
         enddo
      endif

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end
