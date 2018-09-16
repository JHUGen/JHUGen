      subroutine qq_WWqq(p,msq)
      implicit none
      include 'types.f'

c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W(p3,p4)+W(p5,p6)+q(p7)+q(p8);
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
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
      integer:: h1,h2
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,mult
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2),cdotpr,
     & k7341(4),k8341(4),k8342(4),
     & s7341,s8341,s8342
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
     & j8_34_2z(2,4),j8_34_2g(2,4),j7_56_1z(2,4),j7_56_1g(2,4)
      logical:: doHO,doBO
      parameter(spinavge=0.25d0,stat=0.5d0,nfinc=4)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO,mult
      include 'cplx.h'

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
       mult=1d0
c--- rescaling factor for Higgs amplitudes, if anomalous Higgs width
       if (anom_Higgs) then
         mult=chi_higgs**2
       endif
       first=.false.
       call flush(6)
      endif

      if (doHO) then
        Hbit=mult*cone
        Bbit=czip
      elseif (doBO) then
        Hbit=czip
        Bbit=cone
      else
        Hbit=mult*cone
        Bbit=cone
      endif

c--- note that this is the special ordering to agree with Madgraph
      i3=3
      i4=6
      i5=5
      i6=4
c--- this is the MCFM ordering in process.DAT
c      i3=3
c      i4=4
c      i5=5
c      i6=6

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:)=czip
      ampa(:,:,:)=czip
      ampb(:,:,:)=czip

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


C-----setup for (dqcq_dqcq)
      do h1=1,2
      do h2=1,2

      amp(dqcq_dqcq,h1,h2)=
     & +jmid17(1,2,h1,h2)
     & +jvbf17(1,2,h1,h2)
     & +jtwo17(1,2,h1,h2)+jtwo28(2,1,h2,h1)
     & +jZWZa17(1,2,h1,h2)+jZWZb17(1,2,h1,h2)
     & +jtwodiags17(1,2,h1,h2)

      temp(1,4)=temp(1,4)+esq**6*spinavge
     &   *real(amp(dqcq_dqcq,h1,h2)
     & *conjg(amp(dqcq_dqcq,h1,h2)))

      enddo
      enddo

C-----setup for (dqcq_uqsq)
      amp(dqcq_uqsq,1,1)=
     & +jtwoWexch17(1,2)+jtwoWexch28(2,1)
     & +jmidWW17
     & +cdotpr(j7_34_1g(1,:),j8_56_2g(2,:))/s7341
     & +(cdotpr(j7_34_1z(1,:),j8_56_2z(2,:))
     &  -cdotpr(j7_34_1z(1,:),k7341(:))
     &  *cdotpr(k7341(:),j8_56_2z(2,:))/czmass2)/(s7341-czmass2)

      tempw(1,4)=tempw(1,4)+esq**6*spinavge
     &   *real(amp(dqcq_uqsq,1,1)
     & *conjg(amp(dqcq_uqsq,1,1)))

C-----setup for (uqcq_uqcq)
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(uqcq_uqcq,h1,h2)=
     & +jmid17(2,2,h1,h2)
     & +jvbf17(2,2,h1,h2)
     & +jtwo17(2,2,h1,h2)+jtwo28(2,2,h2,h1)
     & +jZWZa17(2,2,h1,h2)+jZWZb17(2,2,h1,h2)
     & +jtwodiags17(2,2,h1,h2)

      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *real(amp(uqcq_uqcq,h1,h2)
     & *conjg(amp(uqcq_uqcq,h1,h2)))

      enddo
      enddo


C-----setup for (dqsq_dqsq)
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(dqsq_dqsq,h1,h2)=
     & +jmid17(1,1,h1,h2)
     & +jvbf17(1,1,h1,h2)
     & +jtwo17(1,1,h1,h2)+jtwo28(1,1,h2,h1)
     & +jZWZa17(1,1,h1,h2)+jZWZb17(1,1,h1,h2)
     & +jtwodiags17(1,1,h1,h2)

      temp(1,3)=temp(1,3)+esq**6*spinavge
     &   *real(amp(dqsq_dqsq,h1,h2)
     & *conjg(amp(dqsq_dqsq,h1,h2)))

      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

C-----setup for (dqdq_dqdq)
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
c-------- ampa
      ampa(dqdq_dqdq,h1,h2)=amp(dqsq_dqsq,h1,h2)

c-------- ampb
      ampb(dqdq_dqdq,h1,h2)=
     & +jmid18(1,1,h1,h2)
     & +jvbf18(1,1,h1,h2)
     & +jtwo18(1,1,h1,h2)+jtwo27(1,1,h2,h1)
     & +jZWZa18(1,1,h1,h2)+jZWZb18(1,1,h1,h2)
     & +jtwodiags18(1,1,h1,h2)

      temp(1,1)=temp(1,1)+esq**6*spinavge
     &   *real(ampa(dqdq_dqdq,h1,h2)
     & *conjg(ampa(dqdq_dqdq,h1,h2)))
      temp(1,1)=temp(1,1)+esq**6*spinavge
     &   *real(ampb(dqdq_dqdq,h1,h2)
     & *conjg(ampb(dqdq_dqdq,h1,h2)))
      if (h1 == h2) then
      temp(1,1)=temp(1,1)-2d0/xn*esq**6*spinavge
     &   *real(ampa(dqdq_dqdq,h1,h2)
     & *conjg(ampb(dqdq_dqdq,h1,h2)))
      endif

      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)

C-----setup for (uquq_uquq)
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
c-------- ampa
      ampa(uquq_uquq,h1,h2)=amp(uqcq_uqcq,h1,h2)

c-------- ampb
      ampb(uquq_uquq,h1,h2)=
     & +jmid18(2,2,h1,h2)
     & +jvbf18(2,2,h1,h2)
     & +jtwo18(2,2,h1,h2)+jtwo27(2,2,h2,h1)
     & +jZWZa18(2,2,h1,h2)+jZWZb18(2,2,h1,h2)
     & +jtwodiags18(2,2,h1,h2)

      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampa(uquq_uquq,h1,h2)
     & *conjg(ampa(uquq_uquq,h1,h2)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampb(uquq_uquq,h1,h2)
     & *conjg(ampb(uquq_uquq,h1,h2)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)-2d0/xn*esq**6*spinavge
     &   *real(ampa(uquq_uquq,h1,h2)
     & *conjg(ampb(uquq_uquq,h1,h2)))
      endif

      enddo
      enddo
      temp(4,4)=temp(2,2)

C-----setup for (dquq_dquq)
c-------- ampb
      ampb(dquq_dquq,1,1)=
     & +jtwoWexch18(1,2)+jtwoWexch27(2,1)
     & +jmidWW18
     & +cdotpr(j8_34_1g(1,:),j7_56_2g(2,:))/s8341
     & +(cdotpr(j8_34_1z(1,:),j7_56_2z(2,:))
     &  -cdotpr(j8_34_1z(1,:),k8341(:))
     &  *cdotpr(k8341(:),j7_56_2z(2,:))/czmass2)/(s8341-czmass2)

      do h1=1,2
      do h2=1,2
c-------- ampa
      ampa(dquq_dquq,h1,h2)=amp(dqcq_dqcq,h1,h2)

      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *real(ampa(dquq_dquq,h1,h2)
     & *conjg(ampa(dquq_dquq,h1,h2)))
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *real(ampb(dquq_dquq,h1,h2)
     & *conjg(ampb(dquq_dquq,h1,h2)))
      if (h1 == h2) then
      temp(1,2)=temp(1,2)-2d0/xn*esq**6*spinavge
     &   *real(ampa(dquq_dquq,h1,h2)
     & *conjg(ampb(dquq_dquq,h1,h2)))
      endif

      enddo
      enddo
      temp(3,4)=temp(1,2)

C-----setup for (uqbq_uqbq)
      do h1=1,2
      do h2=1,2

      amp(uqbq_uqbq,h1,h2)=
     & +jmid17(2,1,h1,h2)
     & +jvbf17(2,1,h1,h2)
     & +jtwo17(2,1,h1,h2)+jtwo28(1,2,h2,h1)
     & +jZWZa17(2,1,h1,h2)+jZWZb17(2,1,h1,h2)
     & +jtwodiags17(2,1,h1,h2)

      temp(2,5)=temp(2,5)+esq**6*spinavge
     &   *real(amp(uqbq_uqbq,h1,h2)
     & *conjg(amp(uqbq_uqbq,h1,h2)))

      enddo
      enddo
      temp(4,5)=temp(2,5)
      temp(2,3)=temp(2,5)

C-----setup for (uqsq_dqcq)
      amp(uqsq_dqcq,1,1)=
     & +jtwoWexch17(2,1)+jtwoWexch28(1,2)
     & +jmidWW28
     & +cdotpr(j8_34_2g(1,:),j7_56_1g(2,:))/s8342
     & +(cdotpr(j8_34_2z(1,:),j7_56_1z(2,:))
     &  -cdotpr(j8_34_2z(1,:),k8342(:))
     &  *cdotpr(k8342(:),j7_56_1z(2,:))/czmass2)/(s8342-czmass2)

      tempw(2,3)=tempw(2,3)+esq**6*spinavge
     &   *real(amp(uqsq_dqcq,1,1)
     & *conjg(amp(uqsq_dqcq,1,1)))

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

