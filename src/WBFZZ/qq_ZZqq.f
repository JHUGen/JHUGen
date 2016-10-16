      subroutine qq_ZZqq_perpart(p,pid,msq)
      implicit none
c--- Author: R.K. Ellis, October 2014
c--- q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
!      include 'first.f'
      include 'spinzerohiggs_anomcoupl.f'
      include 'WWbits.f'
      integer nmax,jmax
      parameter(jmax=12,nmax=10)
      integer j,k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2,h3,h5
      integer pid(mxpart)
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,mult
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),cdotpr,
     & propw71,propw81,propw72,propw82,
     & propw7341,propw7561,propw7342,propw7562,
     & ZZ7341(2,2,2,2,2,2),ZZ7561(2,2,2,2,2,2),
     & ZZ8341(2,2,2,2,2,2),ZZ8561(2,2,2,2,2,2),
     & WWp7341(2),WWm7341(2),WWp7561(2),WWm7561(2),
     & WWp8341(2),WWm8341(2),WWp8561(2),WWm8561(2),
     & j7_1(4,2),j7_2(4,2),j8_1(4,2),j8_2(4,2),
     & j7_34_1(4,2,2,2),j7_34_2(4,2,2,2),
     & j7_56_1(4,2,2,2),j7_56_2(4,2,2,2),
     & j8_34_1(4,2,2,2),j8_34_2(4,2,2,2),
     & j8_56_1(4,2,2,2),j8_56_2(4,2,2,2),
     & jl7_34_1(4,2,2,2),jl7_34_2(4,2,2,2),
     & jl7_56_1(4,2,2,2),jl7_56_2(4,2,2,2),
     & jl8_34_1(4,2,2,2),jl8_34_2(4,2,2,2),
     & jl8_56_1(4,2,2,2),jl8_56_2(4,2,2,2),
     & jw7_34_1(4,2,2),jw7_34_2(4,2,2),
     & jw7_56_1(4,2,2),jw7_56_2(4,2,2),
     & jw8_34_1(4,2,2),jw8_34_2(4,2,2),
     & jw8_56_1(4,2,2),jw8_56_2(4,2,2),
     & j7_3456_1(4,2,2,2,2),j8_3456_2(4,2,2,2,2),
     & j7_3456_2(4,2,2,2,2),j8_3456_1(4,2,2,2,2),
     & jw7_3456_1(4,2,2,2),jw8_3456_2(4,2,2,2),
     & jw7_3456_2(4,2,2,2),jw8_3456_1(4,2,2,2),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),
     & ZZHamp71_82(2,2,2,2,2,2),ZZHamp81_72(2,2,2,2,2,2),
     & WWZZ71_82amp(2,2),WWZZ81_72amp(2,2),
     & WWZZ82_71amp(2,2),WWZZ72_81amp(2,2),
     & gmZ7341(2,2,2,2),gmZ7561(2,2,2,2),gmZ71(2,2,2,2),gmZ82(2,2,2,2),
     & gmZ7342(2,2,2,2),gmZ7562(2,2,2,2),gmZ72(2,2,2,2),gmZ81(2,2,2,2),
     & ll7341(2,2),ll7561(2,2),ll7342(2,2),ll7562(2,2),
     & gmZl7341(2,2,2),gmZl7561(2,2,2),gmZl7342(2,2,2),gmZl7562(2,2,2),
     & gmZl8562(2,2,2),gmZl8342(2,2,2),gmZl8561(2,2,2),gmZl8341(2,2,2),
     & k7341(4),k1567(4),k7342(4),k7562(4),ggWW(2,2),
     & srWWZZ71_82amp(2,2),srWWZZ81_72amp(2,2),
     & srWWZZ82_71amp(2,2),srWWZZ72_81amp(2,2)
c     ,j3_4(4,2),j5_6(4,2),
      logical doHO,doBO
      parameter(spinavge=0.25d0,stat=0.5d0)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO,mult



!$omp threadprivate(doHO,doBO,mult)
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


C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

C--   MARKUS: adding switches to remove VH or VBF contributions
      if( (vvhvvtoggle_vbfvh.eq.0) .and. (j.ge.9) ) cycle ! No VH-like diagram
      if( (vvhvvtoggle_vbfvh.eq.1) .and. (j.le.8) ) cycle ! No VBF-like diagram


c--- propagators and currents are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then
      call setupzprops(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & gmZ7341,gmZ7561,gmZ71,gmZ82,
     & gmZ7342,gmZ7562,gmZ72,gmZ81,
     & ggWW,propw71,propw81,propw72,propw82,
     & propw7341,propw7561,propw7342,propw7562,
     & ll7341,ll7561,ll7342,ll7562,
     & gmZl7341,gmZl7561,gmZl7342,gmZl7562,
     & gmZl8562,gmZl8342,gmZl8561,gmZl8341)

      call jzero(j7(j),j1(j),zab,zba,j7_1)
      call jzero(j7(j),j2(j),zab,zba,j7_2)
      call jzero(j8(j),j1(j),zab,zba,j8_1)
      call jzero(j8(j),j2(j),zab,zba,j8_2)

      call jone(j7(j),3,4,j1(j),za,zb,zab,zba,j7_34_1,jw7_34_1,jl7_34_1)
      call jone(j7(j),3,4,j2(j),za,zb,zab,zba,j7_34_2,jw7_34_2,jl7_34_2)
      call jone(j7(j),5,6,j1(j),za,zb,zab,zba,j7_56_1,jw7_56_1,jl7_56_1)
      call jone(j7(j),5,6,j2(j),za,zb,zab,zba,j7_56_2,jw7_56_2,jl7_56_2)
      call jone(j8(j),3,4,j1(j),za,zb,zab,zba,j8_34_1,jw8_34_1,jl8_34_1)
      call jone(j8(j),3,4,j2(j),za,zb,zab,zba,j8_34_2,jw8_34_2,jl8_34_2)
      call jone(j8(j),5,6,j1(j),za,zb,zab,zba,j8_56_1,jw8_56_1,jl8_56_1)
      call jone(j8(j),5,6,j2(j),za,zb,zab,zba,j8_56_2,jw8_56_2,jl8_56_2)

      call jtwo(j7(j),3,4,5,6,j1(j),za,zb,zab,zba,j7_3456_1,jw7_3456_1)
      call jtwo(j7(j),3,4,5,6,j2(j),za,zb,zab,zba,j7_3456_2,jw7_3456_2)
      call jtwo(j8(j),3,4,5,6,j1(j),za,zb,zab,zba,j8_3456_1,jw8_3456_1)
      call jtwo(j8(j),3,4,5,6,j2(j),za,zb,zab,zba,j8_3456_2,jw8_3456_2)

      k7341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k1567(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      k7342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7562(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))

C-----Singly resonant production in VBF style diagrams
      call ZZSingleres(j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,
     & ZZ7341,WWp7341,WWm7341)
      call ZZSingleres(j1(j),j2(j),5,6,3,4,j7(j),j8(j),za,zb,
     & ZZ7561,WWp7561,WWm7561)
      call ZZSingleres(j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,
     & ZZ8341,WWp8341,WWm8341)
      call ZZSingleres(j1(j),j2(j),5,6,3,4,j8(j),j7(j),za,zb,
     & ZZ8561,WWp8561,WWm8561)
      endif



C----ZZ->ZZ scattering with the exchange of a H
      call ZZHZZamp(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & za,zb,ZZHamp71_82)
      call ZZHZZamp(j1(j),j2(j),3,4,5,6,j8(j),j7(j),
     & za,zb,ZZHamp81_72)
C----Four boson vertex + WW->Higgs diagram
      call WWZZ(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & za,zb,WWZZ71_82amp,srWWZZ71_82amp)
      call WWZZ(j1(j),j2(j),3,4,5,6,j8(j),j7(j),
     & za,zb,WWZZ81_72amp,srWWZZ81_72amp)
      call WWZZ(j2(j),j1(j),3,4,5,6,j8(j),j7(j),
     & za,zb,WWZZ82_71amp,srWWZZ82_71amp)
      call WWZZ(j2(j),j1(j),3,4,5,6,j7(j),j8(j),
     & za,zb,WWZZ72_81amp,srWWZZ72_81amp)

C-----setup for (uqbq_uqbq) (2,5)->(2,5)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqbq_uqbq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))*gmZ7341(2,1,h1,h2)
     & +cdotpr(jl7_34_1(:,2,h1,h3),jl8_56_2(:,1,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,2,h1,h3),jl8_56_2(:,1,h2,h5))*gmZl8562(2,h1,h5)
     & +cdotpr(jl7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))*gmZl7341(1,h2,h3)

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))*gmZ7561(2,1,h1,h2)
     & +cdotpr(jl7_56_1(:,2,h1,h5),jl8_34_2(:,1,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),jl8_34_2(:,1,h2,h3))*gmZl8342(2,h1,h3)
     & +cdotpr(jl7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))*gmZl7561(1,h2,h5)

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))*gmZ82(2,1,h1,h2)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,1,h2,h3,h5))*gmZ71(2,1,h1,h2)

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +ZZ7341(2,1,h1,h2,h3,h5)+ZZ7561(2,1,h1,h2,h5,h3)
      endif

      amp(uqbq_uqbq,h1,h2,h3,h5)=Bbit*amp(uqbq_uqbq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,1,h1,h2,h3,h5)

      temp(2,5)=temp(2,5)+esq**6*spinavge
     &   *dble(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *dconjg(amp(uqbq_uqbq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(4,5)=temp(2,5)
C--------------------------------------------------------------------------
C-----setup for (uqcq_uqcq) (2,4)->(2,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      if (doHO .eqv. .false.) then
      amp(uqcq_uqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))*gmZ7341(2,2,h1,h2)
     & +cdotpr(jl7_34_1(:,2,h1,h3),jl8_56_2(:,2,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,2,h1,h3),jl8_56_2(:,2,h2,h5))*gmZl8562(2,h1,h5)
     & +cdotpr(jl7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))*gmZl7341(2,h2,h3)

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))*gmZ7561(2,2,h1,h2)
     & +cdotpr(jl7_56_1(:,2,h1,h5),jl8_34_2(:,2,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),jl8_34_2(:,2,h2,h3))*gmZl8342(2,h1,h3)
     & +cdotpr(jl7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))*gmZl7561(2,h2,h5)

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))*gmZ82(2,2,h1,h2)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))*gmZ71(2,2,h1,h2)

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +ZZ7341(2,2,h1,h2,h3,h5)+ZZ7561(2,2,h1,h2,h5,h3)
      endif

      amp(uqcq_uqcq,h1,h2,h3,h5)=Bbit*amp(uqcq_uqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,2,h1,h2,h3,h5)

      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqcq_uqcq,h1,h2,h3,h5)))

c      if (j .eq. 1) then
c      write(6,*) h1,h2,h3,h5,esq**6*spinavge
c     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
c     & *dconjg(amp(uqcq_uqcq,h1,h2,h3,h5)))
c      endif

      enddo
      enddo
      enddo
      enddo


C-----setup for uqsq_dqcq W diagrams (2,3)->(1,4)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqsq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(jw7_3456_1(:,2,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,1,h3,h5))*0.5d0/propw71/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,2,h3),jw8_56_2(:,1,h5))
     &  -cdotpr(jw7_34_1(:,2,h3),k7341(:))
     &  *cdotpr(k7341(:),jw8_56_2(:,1,h5))/cwmass2)/propw7341

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,2,h5),jw8_34_2(:,1,h3))
     &  -cdotpr(jw7_56_1(:,2,h5),k1567(:))
     &  *cdotpr(k1567(:),jw8_34_2(:,1,h3))/cwmass2)/propw7561

      if (h3 .eq. 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7341(h5)
      if (h5 .eq. 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7561(h3)

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +srWWZZ71_82amp(h3,h5)
      endif

      amp(uqsq_dqcq,h1,h2,h3,h5)=Bbit*amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWZZ82_71amp(h3,h5)

      tempw(2,3)=tempw(2,3)+esq**6*spinavge
     &   *dble(amp(uqsq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqsq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(2,3)=temp(2,5)

C-----setup for dqcq_uqsq (1,4)-->(2,3)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_uqsq,h1,h2,h3,h5)=
     & +cdotpr(jw7_3456_1(:,1,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,2,h3,h5))*0.5d0/propw71/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,1,h3),jw8_56_2(:,2,h5))
     &  -cdotpr(jw7_34_1(:,1,h3),k7341(:))
     &  *cdotpr(k7341(:),jw8_56_2(:,2,h5))/cwmass2)/propw7341

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,1,h5),jw8_34_2(:,2,h3))
     &  -cdotpr(jw7_56_1(:,1,h5),k1567(:))
     &  *cdotpr(k1567(:),jw8_34_2(:,2,h3))/cwmass2)/propw7561

      if (h3 .eq. 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7341(h5)
      if (h5 .eq. 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7561(h3)

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & -srWWZZ71_82amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif

      amp(dqcq_uqsq,h1,h2,h3,h5)=Bbit*amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWZZ71_82amp(h3,h5)

      tempw(1,4)=tempw(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_uqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_uqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

C-----setup for (dqcq_dqcq) (1,4)-->(1,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))*gmZ7341(1,2,h1,h2)
     & +cdotpr(jl7_34_1(:,1,h1,h3),jl8_56_2(:,2,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,1,h1,h3),jl8_56_2(:,2,h2,h5))*gmZl8562(1,h1,h5)
     & +cdotpr(jl7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))*gmZl7341(2,h2,h3)
      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))*gmZ7561(1,2,h1,h2)
     & +cdotpr(jl7_56_1(:,1,h1,h5),jl8_34_2(:,2,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),jl8_34_2(:,2,h2,h3))*gmZl8342(1,h1,h3)
     & +cdotpr(jl7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))*gmZl7561(2,h2,h5)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))*gmZ82(1,2,h1,h2)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))*gmZ71(1,2,h1,h2)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +ZZ7341(1,2,h1,h2,h3,h5)+ZZ7561(1,2,h1,h2,h5,h3)
      endif

      amp(dqcq_dqcq,h1,h2,h3,h5)=Bbit*amp(dqcq_dqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,2,h1,h2,h3,h5)

      temp(1,4)=temp(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

C--------------------------------------------------------
C-----setup for dquq_dquq W diagrams (1,2)-->(1,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if ((h1.eq.1) .and. (h2.eq.1)) then
      if (doHO .eqv. .false.) then
      ampa(dquq_dquq,h1,h2,h3,h5)=
     & +cdotpr(jw8_3456_1(:,1,h3,h5),j7_2(:,h2))*0.5d0/propw72/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +cdotpr(j8_1(:,h1),jw7_3456_2(:,2,h3,h5))*0.5d0/propw81/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_34_1(:,1,h3),jw7_56_2(:,2,h5))
     &  -cdotpr(jw8_34_1(:,1,h3),k7562(:))
     &  *cdotpr(k7562(:),jw7_56_2(:,2,h5))/cwmass2)/propw7562

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_56_1(:,1,h5),jw7_34_2(:,2,h3))
     &  -cdotpr(jw8_56_1(:,1,h5),k7342(:))
     &  *cdotpr(k7342(:),jw7_34_2(:,2,h3))/cwmass2)/propw7342

      if (h3 .eq. 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8341(h5)
      if (h5 .eq. 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8561(h3)

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & -srWWZZ81_72amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif

      ampa(dquq_dquq,h1,h2,h3,h5)=Bbit*ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWZZ81_72amp(h3,h5)
      endif

C--Fill Z exchange diagrams
      ampb(dquq_dquq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)

      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampa(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampb(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)-2d0/xn*esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(3,4)=temp(1,2)
C-----------------------------------------------------------------

C-----setup for (dqsq_dqsq) (1,3)-->(1,3)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqsq_dqsq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))*gmZ7341(1,1,h1,h2)
     & +cdotpr(jl7_34_1(:,1,h1,h3),jl8_56_2(:,1,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,1,h1,h3),jl8_56_2(:,1,h2,h5))*gmZl8562(1,h1,h5)
     & +cdotpr(jl7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))*gmZl7341(1,h2,h3)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))*gmZ7561(1,1,h1,h2)
     & +cdotpr(jl7_56_1(:,1,h1,h5),jl8_34_2(:,1,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),jl8_34_2(:,1,h2,h3))*gmZl8342(1,h1,h3)
     & +cdotpr(jl7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))*gmZl7561(1,h2,h5)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))*gmZ82(1,1,h2,h1)
     & +cdotpr(j8_3456_2(:,1,h2,h3,h5),j7_1(:,h1))*gmZ71(1,1,h1,h2)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +ZZ7341(1,1,h1,h2,h3,h5)+ZZ7561(1,1,h1,h2,h5,h3)
      endif

      amp(dqsq_dqsq,h1,h2,h3,h5)=Bbit*amp(dqsq_dqsq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,1,h1,h2,h3,h5)

      temp(1,3)=temp(1,3)+esq**6*spinavge
     &   *dble(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

      if ((j.eq.2).or.(j.eq.4).or.(j.eq.6).or.(j.eq.8)
     & .or.(j.eq.10).or.(j.eq.12)) go to 100
C-----setup for ((uquq_uquq)  (2,2)-->(2,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C-----------------ampa
      ampa(uquq_uquq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)

C-----------------ampb
      if (doHO .eqv. .false.) then
      ampb(uquq_uquq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))*gmZ7562(2,2,h1,h2)
     & +cdotpr(jl8_34_1(:,2,h1,h3),jl7_56_2(:,2,h2,h5))*ll7562(h3,h5)
     & +cdotpr(j8_34_1(:,2,h1,h3),jl7_56_2(:,2,h2,h5))*gmZl7562(2,h1,h5)
     & +cdotpr(jl8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))*gmZl8341(2,h2,h3)
      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))*gmZ7342(2,2,h1,h2)
     & +cdotpr(jl8_56_1(:,2,h1,h5),jl7_34_2(:,2,h2,h3))*ll7342(h3,h5)
     & +cdotpr(j8_56_1(:,2,h1,h5),jl7_34_2(:,2,h2,h3))*gmZl7342(2,h1,h3)
     & +cdotpr(jl8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))*gmZl8561(2,h2,h5)

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,2,h1,h3,h5),j7_2(:,h2))*gmZ72(2,2,h1,h2)
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,2,h2,h3,h5))*gmZ81(2,2,h1,h2)

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +ZZ8341(2,2,h1,h2,h3,h5)+ZZ8561(2,2,h1,h2,h5,h3)
      endif

      ampb(uquq_uquq,h1,h2,h3,h5)=Bbit*ampb(uquq_uquq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(2,2,h1,h2,h3,h5)

      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampa(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *dble(ampb(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*esq**6*spinavge
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(4,4)=temp(2,2)


C-----setup for ((dqdq_dqdq)  (1,1)-->(1,1)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

C-----------------ampa
      ampa(dqdq_dqdq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)

C-----------------ampb
      if (doHO .eqv. .false.) then
      ampb(dqdq_dqdq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))*gmZ7562(1,1,h1,h2)
     & +cdotpr(jl8_34_1(:,1,h1,h3),jl7_56_2(:,1,h2,h5))*ll7562(h3,h5)
     & +cdotpr(j8_34_1(:,1,h1,h3),jl7_56_2(:,1,h2,h5))*gmZl7562(1,h1,h5)
     & +cdotpr(jl8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))*gmZl8341(1,h2,h3)
      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))*gmZ7342(1,1,h1,h2)
     & +cdotpr(jl8_56_1(:,1,h1,h5),jl7_34_2(:,1,h2,h3))*ll7342(h3,h5)
     & +cdotpr(j8_56_1(:,1,h1,h5),jl7_34_2(:,1,h2,h3))*gmZl7342(1,h1,h3)
     & +cdotpr(jl8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))*gmZl8561(1,h2,h5)

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,1,h1,h3,h5),j7_2(:,h2))*gmZ72(1,1,h1,h2)
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,1,h2,h3,h5))*gmZ81(1,1,h1,h2)

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +ZZ8341(1,1,h1,h2,h3,h5)+ZZ8561(1,1,h1,h2,h5,h3)
      endif

      ampb(dqdq_dqdq,h1,h2,h3,h5)=Bbit*ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(1,1,h1,h2,h3,h5)

      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampa(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *dble(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*esq**6*spinavge
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)
 100  continue

      if (j.eq.1) then
      do k=1,nf
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nf
      msq(k,l)=temp(k,l)
      enddo
      enddo
      msq(2,3)=msq(2,3)+tempw(2,3)
      msq(1,4)=msq(1,4)+tempw(1,4)

      elseif (j.eq.2) then
      do k=1,nf
      do l=k+1,nf
      msq(l,k)=temp(k,l)
      enddo
      enddo
      msq(3,2)=msq(3,2)+tempw(2,3)
      msq(4,1)=msq(4,1)+tempw(1,4)

      elseif (j.eq.3) then
      do k=-nf,-1
      msq(k,k)=temp(-k,-k)*stat
      do l=k+1,-1
      msq(k,l)=temp(-l,-k)
      enddo
      enddo
      msq(-3,-2)=msq(-3,-2)+tempw(1,4)
      msq(-4,-1)=msq(-4,-1)+tempw(2,3)

      elseif (j.eq.4) then
      do k=-nf,-1
      do l=k+1,-1
      msq(l,k)=temp(-l,-k)
      enddo
      enddo
      msq(-2,-3)=msq(-2,-3)+tempw(1,4)
      msq(-1,-4)=msq(-1,-4)+tempw(2,3)

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
      msq(-1,3)=msq(-1,3)+tempw(2,3)
      msq(-2,4)=msq(-2,4)+tempw(1,4)

c--- qbar-q
      elseif (j.eq.6) then
      do k=-nf,-1
      do l=1,nf
      if (abs(k) .gt. abs(l)) then
      msq(k,l)=temp(l,-k)
      endif
      enddo
      enddo
      msq(-3,1)=msq(-3,1)+tempw(1,4)
      msq(-4,2)=msq(-4,2)+tempw(2,3)

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
      msq(3,-1)=msq(3,-1)+tempw(2,3)
      msq(4,-2)=msq(4,-2)+tempw(1,4)

c--- q-qbar
      elseif (j.eq.8) then
      do k=-nf,-1
      do l=-nf,-1
      if (abs(k) .lt. abs(l)) then
      msq(-k,l)=temp(-k,-l)
      endif
      enddo
      enddo
      msq(1,-3)=msq(1,-3)+tempw(1,4)
      msq(2,-4)=msq(2,-4)+tempw(2,3)

c--- q-qbar extra pieces
      elseif (j.eq.9) then
      do k=1,nf
      do l=1,nf
      if (k .lt. l) then
      msq(k,-k)=msq(k,-k)+temp(k,l)
      endif
      enddo
      enddo
      msq(1,-2)=msq(1,-2)+tempw(1,4) ! d u~ -> c~ s
      msq(3,-4)=msq(1,-2)
      msq(2,-1)=msq(2,-1)+tempw(2,3)
      msq(4,-3)=msq(2,-1)

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
      msq(-2,1)=msq(-2,1)+tempw(1,4) ! u~ d -> c~ s
      msq(-4,3)=msq(-2,1)
      msq(-1,2)=msq(-1,2)+tempw(2,3) ! d~ u -> s~ c
      msq(-3,4)=msq(-1,2)

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

      subroutine qq_ZZqq(p,msq)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),msq(fn:nf,fn:nf)
      integer pid(mxpart)

      call qq_ZZqq_perpart(p,pid,msq)

      end

