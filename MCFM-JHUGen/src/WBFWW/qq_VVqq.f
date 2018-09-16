      subroutine qq_VVqq(p,msq)
      implicit none
      include 'types.f'
c--- Author: J.M.Campbell, December 2014
c--- q(-p1) + q(-p2) -> e-(p3) e-(p4) nu_e(p5) nu_ebar(p6);
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'first.f'
      include 'WWbits.f'
      integer:: nmax,jmax
      parameter(jmax=12,nmax=10)
      integer:: j,k,l,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq,
     & i3,i4,i5,i6,nfinc
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer:: h1,h2,h3,h5
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & propw71,propw81,propw72,propw82,cdotpr,
     & propw7341,propw7561,propw7342,propw7562,
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),
     & k7341(4),k8341(4),k8342(4),
     & kw7341(4),kw1567(4),kw7342(4),kw7562(4),s7341,s8341,s8342
      complex(dp)::
     & jmid17(2,2,2,2),jvbf17(2,2,2,2),jtwodiags17(2,2,2,2),
     &jtwo17(2,2,2,2),jtwo28(2,2,2,2),jZWZa17(2,2,2,2),jZWZb17(2,2,2,2),
     & jmid18(2,2,2,2),jvbf18(2,2,2,2),jtwodiags18(2,2,2,2),
     &jtwo18(2,2,2,2),jtwo27(2,2,2,2),jZWZa18(2,2,2,2),jZWZb18(2,2,2,2),
     & jtwoWexch17(2,2),jtwoWexch28(2,2),
     & j7_34_1z(2,4),j7_34_1g(2,4),j8_56_2z(2,4),j8_56_2g(2,4),
     & jtwoWexch18(2,2),jtwoWexch27(2,2),
     & j8_34_1z(2,4),j8_34_1g(2,4),j7_56_2z(2,4),j7_56_2g(2,4),
     & jmidWW17,jmidWW18,jmidWW28,
     & j8_34_2z(2,4),j8_34_2g(2,4),j7_56_1z(2,4),j7_56_1g(2,4),
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
     & ZZHamp71_82(2,2,2,2,2,2),ZZHamp81_72(2,2,2,2,2,2),
     & WWZZ71_82amp(2,2),WWZZ81_72amp(2,2),
     & gmZ7341(2,2,2,2),gmZ7561(2,2,2,2),gmZ71(2,2,2,2),gmZ82(2,2,2,2),
     & gmZ7342(2,2,2,2),gmZ7562(2,2,2,2),gmZ72(2,2,2,2),gmZ81(2,2,2,2),
     & ll7341(2,2),ll7561(2,2),ll7342(2,2),ll7562(2,2),
     & gmZl7341(2,2,2),gmZl7561(2,2,2),gmZl7342(2,2,2),gmZl7562(2,2,2),
     & gmZl8562(2,2,2),gmZl8342(2,2,2),gmZl8561(2,2,2),gmZl8341(2,2,2),
     & ggWW(2,2),srWWZZ71_82amp(2,2),srWWZZ81_72amp(2,2)
      logical:: doHO,doBO
      parameter(spinavge=0.25d0,stat=0.5d0,nfinc=4)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO

      msq(:,:)=0d0

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      if (first) then
       cwmass2=cplx2(wmass**2,-wmass*wwidth)
       czmass2=cplx2(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
c       cxw=cplx2(xw,0d0) ! DEBUG: Madgraph comparison
       write(6,*)
       write(6,*) '**************** Complex-mass scheme ***************'
       write(6,*) '*                                                  *'
       write(6,77) cwmass2
       write(6,78) czmass2
       write(6,79) cxw
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       write(6,*)
       doHO=.false.
       doBO=.false.
       if     (runstring(4:5) == 'HO') then
         doHO=.true.
       write(6,*) '>>>>>>>>>>>>>> Higgs contribution only <<<<<<<<<<<<<'
       write(6,*)
       elseif (runstring(4:5) == 'BO') then
         doBO=.true.
       write(6,*)
       write(6,*) '>>>>>>>>>>> Background contribution only <<<<<<<<<<<'
       write(6,*)
       endif
       call flush(6)
       first=.false.
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

c--- note that this is the ordering in process.DAT and in Madgraph
      i3=3
      i4=6
      i5=5
      i6=4

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      k7341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j7(j),:,j7(j)))
      k8341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      k8342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      s7341=cdotpr(k7341,k7341)
      s8341=cdotpr(k8341,k8341)
      s8342=cdotpr(k8342,k8342)

c-----------------------------------------
c---- SET-UP FOR WW-like AMPLITUDES ------
c-----------------------------------------

c--- These contributions contain Hbit and Bbit
c--- contribution from jVBF
      call ampvbf(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jvbf17)
      call ampvbf(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jvbf18)

c--- W mid diagrams: contribution from jmidWW
      call ampmidWW(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jmidWW17)
      call ampmidWW(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jmidWW18)
      call ampmidWW(j2(j),j1(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jmidWW28)

c--- these are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then

c--- mid diagrams: contribution from jcentre
      call ampmid(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jmid17)
      call ampmid(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jmid18)
c--- contribution from jtwoWW
      call amp2current(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jtwo17)
      call amp2current(j2(j),j1(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jtwo28)
      call amp2current(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jtwo18)
      call amp2current(j2(j),j1(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jtwo27)

c--- Z/W/Z diagrams: contribution from jZWZ
      call ampZWZ(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jZWZa17)
      call ampZWZ(j1(j),j2(j),i5,i6,i3,i4,j7(j),j8(j),za,zb,jZWZb17)
      call ampZWZ(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jZWZa18)
      call ampZWZ(j1(j),j2(j),i5,i6,i3,i4,j8(j),j7(j),za,zb,jZWZb18)

c--- W-exchange diagrams: contribution from jtwodiags
      call amptwodiags(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,
     & jtwodiags17)
      call amptwodiags(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,
     & jtwodiags18)

c--- W exchange diagrams for flavor-changing contributions: contribution from jtwoWexch
      call amp2currentw(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,
     & jtwoWexch17)
      call amp2currentw(j2(j),j1(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,
     & jtwoWexch28)
      call amp2currentw(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,
     & jtwoWexch18)
      call amp2currentw(j2(j),j1(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,
     & jtwoWexch27)

      call jonew(j7(j),i3,i4,j1(j),za,zb,zab,j7_34_1z,j7_34_1g)
      call jonew(j8(j),i5,i6,j2(j),za,zb,zab,j8_56_2z,j8_56_2g)
      call jonew(j8(j),i3,i4,j1(j),za,zb,zab,j8_34_1z,j8_34_1g)
      call jonew(j7(j),i5,i6,j2(j),za,zb,zab,j7_56_2z,j7_56_2g)
      call jonew(j8(j),i3,i4,j2(j),za,zb,zab,j8_34_2z,j8_34_2g)
      call jonew(j7(j),i5,i6,j1(j),za,zb,zab,j7_56_1z,j7_56_1g)

      else

      jmid17=czip
      jmid18=czip
      jtwo17=czip
      jtwo28=czip
      jtwo18=czip
      jtwo27=czip
      jZWZa17=czip
      jZWZb17=czip
      jZWZa18=czip
      jZWZb18=czip
      jtwodiags17=czip
      jtwodiags18=czip
      jtwoWexch17=czip
      jtwoWexch28=czip
      jtwoWexch18=czip
      jtwoWexch27=czip
      j7_34_1z=czip
      j7_34_1g=czip
      j8_56_2z=czip
      j8_56_2g=czip
      j8_34_1z=czip
      j8_34_1g=czip
      j7_56_2z=czip
      j7_56_2g=czip
      j8_34_2z=czip
      j8_34_2g=czip
      j7_56_1z=czip
      j7_56_1g=czip

      endif

c-----------------------------------------
c---- SET-UP FOR ZZ-like AMPLITUDES ------
c-----------------------------------------

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

      kw7341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      kw1567(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      kw7342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      kw7562(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(5,:,5)
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


C-----setup for (dqcq_dqcq)
c----- WW-like
      do h1=1,2
      do h2=1,2

      amp(dqcq_dqcq,h1,h2,1,1)=
     & +jmid17(1,2,h1,h2)
     & +jvbf17(1,2,h1,h2)
     & +jtwo17(1,2,h1,h2)+jtwo28(2,1,h2,h1)
     & +jZWZa17(1,2,h1,h2)+jZWZb17(1,2,h1,h2)
     & +jtwodiags17(1,2,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(dqcq_dqcq,h1,h2,1,1)=-amp(dqcq_dqcq,h1,h2,1,1)

      enddo
      enddo

c----- ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
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

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,2,h1,h2,h3,h5)

      temp(1,4)=temp(1,4)+esq**6*spinavge
     &   *real(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *conjg(amp(dqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo


C-----setup for (dqcq_uqsq)
c----- WW-like
      amp(dqcq_uqsq,1,1,1,1)=
     & +jtwoWexch17(1,2)+jtwoWexch28(2,1)
     & +jmidWW17
     & +cdotpr(j7_34_1g(1,:),j8_56_2g(2,:))/s7341
     & +(cdotpr(j7_34_1z(1,:),j8_56_2z(2,:))
     &  -cdotpr(j7_34_1z(1,:),k7341(:))
     &  *cdotpr(k7341(:),j8_56_2z(2,:))/czmass2)/(s7341-czmass2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(dqcq_uqsq,1,1,1,1)=-amp(dqcq_uqsq,1,1,1,1)

c----- ZZ-like
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +cdotpr(jw7_3456_1(:,1,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,2,h3,h5))*0.5d0/propw71/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,1,h3),jw8_56_2(:,2,h5))
     &  -cdotpr(jw7_34_1(:,1,h3),kw7341(:))
     &  *cdotpr(kw7341(:),jw8_56_2(:,2,h5))/cwmass2)/propw7341

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,1,h5),jw8_34_2(:,2,h3))
     &  -cdotpr(jw7_56_1(:,1,h5),kw1567(:))
     &  *cdotpr(kw1567(:),jw8_34_2(:,2,h3))/cwmass2)/propw7561

      if (h3 == 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7341(h5)
      if (h5 == 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7561(h3)

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & -srWWZZ71_82amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWZZ71_82amp(h3,h5)

      tempw(1,4)=tempw(1,4)+esq**6*spinavge
     &   *real(amp(dqcq_uqsq,h1,h2,h3,h5)
     & *conjg(amp(dqcq_uqsq,h1,h2,h3,h5)))

      enddo
      enddo
      enddo
      enddo

C-----setup for (uqcq_uqcq)
c----- WW-like
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(uqcq_uqcq,h1,h2,1,1)=
     & +jmid17(2,2,h1,h2)
     & +jvbf17(2,2,h1,h2)
     & +jtwo17(2,2,h1,h2)+jtwo28(2,2,h2,h1)
     & +jZWZa17(2,2,h1,h2)+jZWZb17(2,2,h1,h2)
     & +jtwodiags17(2,2,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(uqcq_uqcq,h1,h2,1,1)=-amp(uqcq_uqcq,h1,h2,1,1)

      enddo
      enddo

c----- ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      if (doHO .eqv. .false.) then
      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
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

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,2,h1,h2,h3,h5)

      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *real(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_uqcq,h1,h2,h3,h5)))

      enddo
      enddo
      enddo
      enddo


C-----setup for (dqsq_dqsq)
c----- WW-like
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(dqsq_dqsq,h1,h2,1,1)=
     & +jmid17(1,1,h1,h2)
     & +jvbf17(1,1,h1,h2)
     & +jtwo17(1,1,h1,h2)+jtwo28(1,1,h2,h1)
     & +jZWZa17(1,1,h1,h2)+jZWZb17(1,1,h1,h2)
     & +jtwodiags17(1,1,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(dqsq_dqsq,h1,h2,1,1)=-amp(dqsq_dqsq,h1,h2,1,1)

      enddo
      enddo

c----- ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
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

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,1,h1,h2,h3,h5)

      temp(1,3)=temp(1,3)+esq**6*spinavge
     &   *real(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(dqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

C-----setup for (dqdq_dqdq)
c----- WW-like
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
c-------- ampa
      ampa(dqdq_dqdq,h1,h2,:,:)=amp(dqsq_dqsq,h1,h2,:,:)

c-------- ampb: WW-like
      ampb(dqdq_dqdq,h1,h2,1,1)=
     & +jmid18(1,1,h1,h2)
     & +jvbf18(1,1,h1,h2)
     & +jtwo18(1,1,h1,h2)+jtwo27(1,1,h2,h1)
     & +jZWZa18(1,1,h1,h2)+jZWZb18(1,1,h1,h2)
     & +jtwodiags18(1,1,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      ampb(dqdq_dqdq,h1,h2,1,1)=-ampb(dqdq_dqdq,h1,h2,1,1)

      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

C-------- ampb: ZZ-like
      if (doHO .eqv. .false.) then
      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
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

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(1,1,h1,h2,h3,h5)

      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *real(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampa(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *real(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(1,1)=temp(1,1)-2d0/xn*esq**6*spinavge
     & *real(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)

C-----setup for (uquq_uquq)
c----- WW-like
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
c-------- ampa
      ampa(uquq_uquq,h1,h2,:,:)=amp(uqcq_uqcq,h1,h2,:,:)

c-------- ampb: WW-like
      ampb(uquq_uquq,h1,h2,1,1)=
     & +jmid18(2,2,h1,h2)
     & +jvbf18(2,2,h1,h2)
     & +jtwo18(2,2,h1,h2)+jtwo27(2,2,h2,h1)
     & +jZWZa18(2,2,h1,h2)+jZWZb18(2,2,h1,h2)
     & +jtwodiags18(2,2,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      ampb(uquq_uquq,h1,h2,1,1)=-ampb(uquq_uquq,h1,h2,1,1)

      enddo
      enddo

c-------- ampb: ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      if (doHO .eqv. .false.) then
      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
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

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(2,2,h1,h2,h3,h5)

      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *real(ampa(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *real(ampb(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)-2d0/xn*esq**6*spinavge
     & *real(ampa(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(4,4)=temp(2,2)

C-----setup for (dquq_dquq)
c----- WW-like
c-------- ampb
      ampa(dquq_dquq,1,1,1,1)=
     & +jtwoWexch18(1,2)+jtwoWexch27(2,1)
     & +jmidWW18
     & +cdotpr(j8_34_1g(1,:),j7_56_2g(2,:))/s8341
     & +(cdotpr(j8_34_1z(1,:),j7_56_2z(2,:))
     &  -cdotpr(j8_34_1z(1,:),k8341(:))
     &  *cdotpr(k8341(:),j7_56_2z(2,:))/czmass2)/(s8341-czmass2)

c--- the following minus sign accounts for one less fermion loop in interference
      ampa(dquq_dquq,1,1,1,1)=-ampa(dquq_dquq,1,1,1,1)


      do h1=1,2
      do h2=1,2
c-------- ampa: WW-like
      ampb(dquq_dquq,h1,h2,:,:)=amp(dqcq_dqcq,h1,h2,:,:)

      enddo
      enddo

c-------- ampa: ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if ((h1==1) .and. (h2==1)) then
      if (doHO .eqv. .false.) then
      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +cdotpr(jw8_3456_1(:,1,h3,h5),j7_2(:,h2))*0.5d0/propw72/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +cdotpr(j8_1(:,h1),jw7_3456_2(:,2,h3,h5))*0.5d0/propw81/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_34_1(:,1,h3),jw7_56_2(:,2,h5))
     &  -cdotpr(jw8_34_1(:,1,h3),kw7562(:))
     &  *cdotpr(kw7562(:),jw7_56_2(:,2,h5))/cwmass2)/propw7562

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_56_1(:,1,h5),jw7_34_2(:,2,h3))
     &  -cdotpr(jw8_56_1(:,1,h5),kw7342(:))
     &  *cdotpr(kw7342(:),jw7_34_2(:,2,h3))/cwmass2)/propw7342

      if (h3 == 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8341(h5)
      if (h5 == 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8561(h3)

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & -srWWZZ81_72amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWZZ81_72amp(h3,h5)
      endif

      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *real(ampa(dquq_dquq,h1,h2,h3,h5)
     & *conjg(ampa(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *real(ampb(dquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)-2d0/xn*esq**6*spinavge
     &   *real(ampa(dquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

      temp(3,4)=temp(1,2)

C-----setup for (uqbq_uqbq)
c----- WW-like
      do h1=1,2
      do h2=1,2

      amp(uqbq_uqbq,h1,h2,1,1)=
     & +jmid17(2,1,h1,h2)
     & +jvbf17(2,1,h1,h2)
     & +jtwo17(2,1,h1,h2)+jtwo28(1,2,h2,h1)
     & +jZWZa17(2,1,h1,h2)+jZWZb17(2,1,h1,h2)
     & +jtwodiags17(2,1,h1,h2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(uqbq_uqbq,h1,h2,1,1)=-amp(uqbq_uqbq,h1,h2,1,1)


      enddo
      enddo

c----- ZZ-like
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
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

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,1,h1,h2,h3,h5)

      temp(2,5)=temp(2,5)+esq**6*spinavge
     &   *real(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *conjg(amp(uqbq_uqbq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(4,5)=temp(2,5)
      temp(2,3)=temp(2,5)

C-----setup for (uqsq_dqcq)
c----- WW-like
      amp(uqsq_dqcq,1,1,1,1)=
     & +jtwoWexch17(2,1)+jtwoWexch28(1,2)
     & +jmidWW28
     & +cdotpr(j8_34_2g(1,:),j7_56_1g(2,:))/s8342
     & +(cdotpr(j8_34_2z(1,:),j7_56_1z(2,:))
     &  -cdotpr(j8_34_2z(1,:),k8342(:))
     &  *cdotpr(k8342(:),j7_56_1z(2,:))/czmass2)/(s8342-czmass2)

c--- the following minus sign accounts for one less fermion loop in interference
      amp(uqsq_dqcq,1,1,1,1)=-amp(uqsq_dqcq,1,1,1,1)

c----- ZZ-like
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +cdotpr(jw7_3456_1(:,2,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,1,h3,h5))*0.5d0/propw71/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,2,h3),jw8_56_2(:,1,h5))
     &  -cdotpr(jw7_34_1(:,2,h3),kw7341(:))
     &  *cdotpr(kw7341(:),jw8_56_2(:,1,h5))/cwmass2)/propw7341

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,2,h5),jw8_34_2(:,1,h3))
     &  -cdotpr(jw7_56_1(:,2,h5),kw1567(:))
     &  *cdotpr(kw1567(:),jw8_34_2(:,1,h3))/cwmass2)/propw7561

      if (h3 == 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7341(h5)
      if (h5 == 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7561(h3)

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +srWWZZ71_82amp(h3,h5)
      endif

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWZZ71_82amp(h3,h5)

      tempw(2,3)=tempw(2,3)+esq**6*spinavge
     &   *real(amp(uqsq_dqcq,h1,h2,h3,h5)
     & *conjg(amp(uqsq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

c--- fill matrix elements
      if (j == 1) then

      do k=1,nfinc
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nfinc
      msq(k,l)=temp(k,l)
      enddo
      enddo
      msq(2,3)=msq(2,3)+tempw(2,3)
      msq(1,4)=msq(1,4)+tempw(1,4)

      elseif (j==2) then
      do k=1,nfinc
      do l=k+1,nfinc
      msq(l,k)=temp(k,l)
      enddo
      enddo
      msq(3,2)=msq(3,2)+tempw(2,3)
      msq(4,1)=msq(4,1)+tempw(1,4)

      elseif (j==3) then
      do k=-nfinc,-1
      msq(k,k)=temp(-k,-k)*stat
      do l=k+1,-1
      msq(k,l)=temp(-l,-k)
      enddo
      enddo
      msq(-3,-2)=msq(-3,-2)+tempw(1,4)
      msq(-4,-1)=msq(-4,-1)+tempw(2,3)

      elseif (j==4) then
      do k=-nfinc,-1
      do l=k+1,-1
      msq(l,k)=temp(-l,-k)
      enddo
      enddo
      msq(-2,-3)=msq(-2,-3)+tempw(1,4)
      msq(-1,-4)=msq(-1,-4)+tempw(2,3)

c--- qbar-q
      elseif (j==5) then
      do k=-nfinc,-1
      msq(k,-k)=temp(-k,-k)
      do l=1,nfinc
      if (abs(k) < abs(l)) then
      msq(k,l)=temp(-k,l)
      endif
      enddo
      enddo
      msq(-1,3)=msq(-1,3)+tempw(2,3)
      msq(-2,4)=msq(-2,4)+tempw(1,4)

c--- qbar-q
      elseif (j==6) then
      do k=-nfinc,-1
      do l=1,nfinc
      if (abs(k) > abs(l)) then
      msq(k,l)=temp(l,-k)
      endif
      enddo
      enddo
      msq(-3,1)=msq(-3,1)+tempw(1,4)
      msq(-4,2)=msq(-4,2)+tempw(2,3)

c--- q-qbar
      elseif (j==7) then
      do k=-nfinc,-1
      msq(-k,k)=temp(-k,-k)
      do l=1,nfinc
      if (abs(k) < abs(l)) then
      msq(l,k)=temp(-k,l)
      endif
      enddo
      enddo
      msq(3,-1)=msq(3,-1)+tempw(2,3)
      msq(4,-2)=msq(4,-2)+tempw(1,4)

c--- q-qbar
      elseif (j==8) then
      do k=-nfinc,-1
      do l=-nfinc,-1
      if (abs(k) < abs(l)) then
      msq(-k,l)=temp(-k,-l)
      endif
      enddo
      enddo
      msq(1,-3)=msq(1,-3)+tempw(1,4)
      msq(2,-4)=msq(2,-4)+tempw(2,3)

c--- q-qbar extra pieces
      elseif (j==9) then
      do k=1,nfinc
      do l=1,nfinc
      if (k < l) then
      msq(k,-k)=msq(k,-k)+temp(k,l)
      endif
      enddo
      enddo
      msq(1,-2)=msq(1,-2)+tempw(1,4) ! d u~ -> c~ s
      msq(3,-4)=msq(1,-2)
      msq(2,-1)=msq(2,-1)+tempw(2,3)
      msq(4,-3)=msq(2,-1)

c--- q-qbar extra pieces
      elseif (j==10) then
      do k=1,nfinc
      do l=1,nfinc
      if (k > l) then
      msq(k,-k)=msq(k,-k)+temp(l,k)
      endif
      enddo
      enddo

c--- qbar-q extra pieces
      elseif (j==11) then
      do k=1,nfinc
      do l=1,nfinc
      if (k < l) then
      msq(-k,k)=msq(-k,k)+temp(k,l)
      endif
      enddo
      enddo
      msq(-2,1)=msq(-2,1)+tempw(1,4) ! u~ d -> c~ s
      msq(-4,3)=msq(-2,1)
      msq(-1,2)=msq(-1,2)+tempw(2,3) ! d~ u -> s~ c
      msq(-3,4)=msq(-1,2)

c--- qbar-q extra pieces
      elseif (j==12) then
      do k=1,nfinc
      do l=1,nfinc
      if (k > l) then
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

