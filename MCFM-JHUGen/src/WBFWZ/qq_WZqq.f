      subroutine qq_WZqq(pin,msq)
      implicit none
      include 'types.f'

c--- Author: J.M. Campbell, January 2015
c--- q(-p1)+q(-p2)->W(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      include 'first.f'
      include 'nwz.f'
      include 'WWbits.f'
      integer:: nmax,jmax
      parameter(jmax=8,nmax=7)
      integer:: j,k,
     & uqsq_dqsq,uqcq_dqcq,uqdq_dqdq,cqdq_dqsq,
     & uqcq_uqsq,uquq_dquq,uquq_uqdq
      parameter(
     & uqsq_dqsq=1,uqcq_dqcq=2,uqdq_dqdq=3,cqdq_dqsq=4,
     & uqcq_uqsq=5,uquq_dquq=6,uquq_uqdq=7)
      integer:: h1,h2,h3,h5
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempb(fn:nf,fn:nf),spinavge,pin(mxpart,4)
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),cdotpr,
     & k7341(4),k7561(4),k7342(4),k7562(4),
     & tmpz1(mxpart,mxpart),tmpz2(mxpart,4,mxpart)
      complex(dp)::
     & jz7_56_1z(4,2,2,2),jz7_56_1g(4,2,2,2),jz7_56_1w(4,2,2,2),
     & jz7_56_2z(4,2,2,2),jz7_56_2g(4,2,2,2),jz7_56_2w(4,2,2,2),
     & jz8_56_1z(4,2,2,2),jz8_56_1g(4,2,2,2),jz8_56_1w(4,2,2,2),
     & jz8_56_2z(4,2,2,2),jz8_56_2g(4,2,2,2),jz8_56_2w(4,2,2,2),
     & jw7_34_1z(2,4),jw7_34_1g(2,4),jw7_34_1w(4,2,2),jw7_34_2w(4,2,2),
     & jw7_34_2z(2,4),jw7_34_2g(2,4),
     & jw8_34_1z(2,4),jw8_34_1g(2,4),jw8_34_1w(4,2,2),
     & jw8_34_2z(2,4),jw8_34_2g(2,4),jw8_34_2w(4,2,2),
     & ampW2_1728(2,2),ampZ2_1728(2,2,2),ampmid_1728(2,2,2),
     & ampW2_1827(2,2),ampZ2_1827(2,2,2),ampmid_1827(2,2,2),
     & ampW2_2817(2,2),ampZ2_2817(2,2,2),ampmid_2817(2,2,2),
     & ampW2_2718(2,2),ampZ2_2718(2,2,2),ampmid_2718(2,2,2),
     & s7341,s7561,s7562,s7342
      logical,save:: doHO,doBO
      parameter(spinavge=0.25_dp)
      integer,parameter::j1(jmax)=(/1,2,8,7,7,2,7,1/)
      integer,parameter::j2(jmax)=(/2,1,7,8,2,7,1,7/)
      integer,parameter::j7(jmax)=(/7,7,1,1,1,8,2,8/)
      integer,parameter::j8(jmax)=(/8,8,2,2,8,1,8,2/)
      real(dp), save:: mult
!$omp threadprivate(doHO,doBO,mult)
      msq(:,:)=0._dp

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      if (first) then
       cwmass2=cplx2(wmass**2,-wmass*wwidth)
       czmass2=cplx2(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
c       cxw=cplx2(xw,0._dp) ! DEBUG: Madgraph comparison
c       write(6,*)
c       write(6,*) '**************** Complex-mass scheme ***************'
c       write(6,*) '*                                                  *'
c       write(6,77) cwmass2
c       write(6,78) czmass2
c       write(6,79) cxw
c       write(6,*) '*                                                  *'
c       write(6,*) '****************************************************'
c       write(6,*)
       doHO=.false.
       doBO=.false.
       if     (runstring(4:5) == 'HO') then
         doHO=.true.
c       write(6,*) '>>>>>>>>>>>>>> Higgs contribution only <<<<<<<<<<<<<'
c       write(6,*)
       elseif (runstring(4:5) == 'BO') then
         doBO=.true.
c       write(6,*)
c       write(6,*) '>>>>>>>>>>> Background contribution only <<<<<<<<<<<'
c       write(6,*)
       endif
       mult=1._dp
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

      if (nwz == +1) then
        p(:,:)=pin(:,:)
      else
c--- W- amplitudes obtained by parity transformation on W+
        p(1:2,:)=pin(1:2,:)
        p(7,:)=pin(8,:)
        p(8,:)=pin(7,:)
        p(3,:)=pin(4,:)
        p(4,:)=pin(3,:)
        p(5,:)=pin(6,:)
        p(6,:)=pin(5,:)
      endif

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)
      if (nwz == -1) then
c--- parity transformation corresponds to complex conjugation
        tmpz1(1:8,1:8)=za(1:8,1:8)
        za(1:8,1:8)=zb(1:8,1:8)
        zb(1:8,1:8)=tmpz1(1:8,1:8)
        tmpz2(1:8,:,1:8)=zab(1:8,:,1:8)
        zab(1:8,:,1:8)=zba(1:8,:,1:8)
        zba(1:8,:,1:8)=tmpz2(1:8,:,1:8)
      endif

      do j=1,jmax
      temp(:,:)=0._dp
      tempb(:,:)=0._dp
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

c--- propagators and currents are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then
c--- jz7_34_1z:  current for Z -> e^- e^+ qq
c--- jz7_34_1g:  current for g -> e^- e^+ qq
c--- jz7_34_1w:  current for W -> e^- e^+ qq
      call jonez
     & (j7(j),5,6,j1(j),za,zb,zab,zba,jz7_56_1z,jz7_56_1g,jz7_56_1w)
      call jonez
     & (j7(j),5,6,j2(j),za,zb,zab,zba,jz7_56_2z,jz7_56_2g,jz7_56_2w)
      call jonez
     & (j8(j),5,6,j1(j),za,zb,zab,zba,jz8_56_1z,jz8_56_1g,jz8_56_1w)
      call jonez
     & (j8(j),5,6,j2(j),za,zb,zab,zba,jz8_56_2z,jz8_56_2g,jz8_56_2w)

c--- jw7_34_1z:  current for Z -> e^- nu~_e qq
c--- jw7_34_1g:  current for g -> e^- nu~_e qq
c---- note: use WZ version that always corresponds to W+ emission
      call jonewWZ(j7(j),3,4,j1(j),za,zb,zab,jw7_34_1z,jw7_34_1g)
      call jonewWZ(j8(j),3,4,j1(j),za,zb,zab,jw8_34_1z,jw8_34_1g)
      call jonewWZ(j8(j),3,4,j2(j),za,zb,zab,jw8_34_2z,jw8_34_2g)
      call jonewWZ(j7(j),3,4,j2(j),za,zb,zab,jw7_34_2z,jw7_34_2g)

c--- jw7_34_1w:  current for W -> e^- nu~_e qq
      call jww(j7(j),3,4,j1(j),za,zb,zab,jw7_34_1w)
      call jww(j7(j),3,4,j2(j),za,zb,zab,jw7_34_2w)
      call jww(j8(j),3,4,j1(j),za,zb,zab,jw8_34_1w)
      call jww(j8(j),3,4,j2(j),za,zb,zab,jw8_34_2w)

c---- W2 contributions
      call jW2exch2(j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,ampW2_1728)
      call jW2exch2(j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,ampW2_1827)
      call jW2exch2(j2(j),j1(j),3,4,5,6,j8(j),j7(j),za,zb,ampW2_2817)
      call jW2exch2(j2(j),j1(j),3,4,5,6,j7(j),j8(j),za,zb,ampW2_2718)

c---- Z2 contributions
      call jZexch2(j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,ampZ2_1728)
      call jZexch2(j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,ampZ2_1827)
      call jZexch2(j2(j),j1(j),3,4,5,6,j8(j),j7(j),za,zb,ampZ2_2817)
      call jZexch2(j2(j),j1(j),3,4,5,6,j7(j),j8(j),za,zb,ampZ2_2718)

      k7341(:)=half*(zab(j1(j),:,j1(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7561(:)=half*(zab(j1(j),:,j1(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      k7342(:)=half*(zab(j2(j),:,j2(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7562(:)=half*(zab(j2(j),:,j2(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      s7341=cdotpr(k7341,k7341)
      s7561=cdotpr(k7561,k7561)
      s7562=cdotpr(k7562,k7562)
      s7342=cdotpr(k7342,k7342)

      else

      ampZ2_1728=czip
      ampZ2_1827=czip
      ampZ2_2817=czip
      ampZ2_2718=czip
      ampW2_1728=czip
      ampW2_1827=czip
      ampW2_2817=czip
      ampW2_2718=czip
      jw7_34_1w=czip
      jw7_34_2w=czip
      jw8_34_1w=czip
      jw8_34_2w=czip
      jw7_34_1z=czip
      jw7_34_1g=czip
      jw8_34_1z=czip
      jw8_34_1g=czip
      jw8_34_2z=czip
      jw8_34_2g=czip
      jw7_34_2z=czip
      jw7_34_2g=czip
      jz7_56_1z=czip
      jz7_56_1g=czip
      jz7_56_1w=czip
      jz7_56_2z=czip
      jz7_56_2g=czip
      jz7_56_2w=czip
      jz8_56_1z=czip
      jz8_56_1g=czip
      jz8_56_1w=czip
      jz8_56_2z=czip
      jz8_56_2g=czip
      jz8_56_2w=czip

      endif

c---- Mid contributions
      call ampmidWZ(j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,ampmid_1728)
      call ampmidWZ(j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,ampmid_1827)
      call ampmidWZ(j2(j),j1(j),3,4,5,6,j8(j),j7(j),za,zb,ampmid_2817)
      call ampmidWZ(j2(j),j1(j),3,4,5,6,j7(j),j8(j),za,zb,ampmid_2718)

C-----setup for (uqsq_dqsq) (2,3)->(1,3)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      amp(uqsq_dqsq,1,h2,1,h5)=
     & +cdotpr(jw7_34_1g(2,:),jz8_56_2g(:,1,h2,h5))/s7341
     & +(cdotpr(jw7_34_1z(2,:),jz8_56_2z(:,1,h2,h5))
     &  -cdotpr(jw7_34_1z(2,:),k7341(:))
     &  *cdotpr(k7341(:),jz8_56_2z(:,1,h2,h5))/czmass2)/(s7341-czmass2)
      amp(uqsq_dqsq,1,h2,1,h5)=amp(uqsq_dqsq,1,h2,1,h5)
     & +(cdotpr(jz7_56_1w(:,2,1,h5),jw8_34_2w(:,1,h2))
     &  -cdotpr(jz7_56_1w(:,2,1,h5),k7561(:))
     &  *cdotpr(k7561(:),jw8_34_2w(:,1,h2))/cwmass2)/(s7561-cwmass2)
      amp(uqsq_dqsq,1,h2,1,h5)=amp(uqsq_dqsq,1,h2,1,h5)
     & +ampZ2_1728(1,h2,h5)
      enddo
      enddo
      do h5=1,2
      amp(uqsq_dqsq,1,1,1,h5)=amp(uqsq_dqsq,1,1,1,h5)
     & +ampW2_1728(1,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      amp(uqsq_dqsq,1,h2,1,h5)=amp(uqsq_dqsq,1,h2,1,h5)
     & +ampmid_1728(1,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,3)=temp(2,3)+esq**6*spinavge
     &   *real(amp(uqsq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(uqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

c-----setup for (cqdq_dqsq)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      amp(cqdq_dqsq,1,h2,1,h5)=
     & +cdotpr(jw8_34_1g(2,:),jz7_56_2g(:,1,h2,h5))/s7562
     & +(cdotpr(jw8_34_1z(2,:),jz7_56_2z(:,1,h2,h5))
     &  -cdotpr(jw8_34_1z(2,:),k7562(:))
     &  *cdotpr(k7562(:),jz7_56_2z(:,1,h2,h5))/czmass2)/(s7562-czmass2)
      amp(cqdq_dqsq,1,h2,1,h5)=amp(cqdq_dqsq,1,h2,1,h5)
     & +(cdotpr(jz8_56_1w(:,2,1,h5),jw7_34_2w(:,1,h2))
     &  -cdotpr(jz8_56_1w(:,2,1,h5),k7342(:))
     &  *cdotpr(k7342(:),jw7_34_2w(:,1,h2))/cwmass2)/(s7342-cwmass2)
      amp(cqdq_dqsq,1,h2,1,h5)=amp(cqdq_dqsq,1,h2,1,h5)
     & +ampZ2_1827(1,h2,h5)
      enddo
      enddo
      do h5=1,2
      amp(cqdq_dqsq,1,1,1,h5)=amp(cqdq_dqsq,1,1,1,h5)
     & +ampW2_1827(1,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      amp(cqdq_dqsq,1,h2,1,h5)=amp(cqdq_dqsq,1,h2,1,h5)
     & +ampmid_1827(1,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(4,1)=temp(4,1)+esq**6*spinavge
     &   *real(amp(cqdq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(cqdq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo


C-----setup for (uqdq_dqdq) (2,1)->(1,1)
      ampa(uqdq_dqdq,:,:,:,:)=amp(uqsq_dqsq,:,:,:,:)
      ampb(uqdq_dqdq,:,:,:,:)=amp(cqdq_dqsq,:,:,:,:)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,1)=temp(2,1)+esq**6*spinavge
     &   *real(ampa(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampa(uqdq_dqdq,h1,h2,h3,h5)))
      temp(2,1)=temp(2,1)+esq**6*spinavge
     &   *real(ampb(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(uqdq_dqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,1)=temp(2,1)-2._dp/xn*esq**6*spinavge
     &   *real(ampa(uqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(uqdq_dqdq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo


C-----setup for (uqcq_dqcq) (2,4)->(1,4)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      amp(uqcq_dqcq,1,h2,1,h5)=
     & +cdotpr(jw7_34_1g(2,:),jz8_56_2g(:,2,h2,h5))/s7341
     & +(cdotpr(jw7_34_1z(2,:),jz8_56_2z(:,2,h2,h5))
     &  -cdotpr(jw7_34_1z(2,:),k7341(:))
     &  *cdotpr(k7341(:),jz8_56_2z(:,2,h2,h5))/czmass2)/(s7341-czmass2)
      amp(uqcq_dqcq,1,h2,1,h5)=amp(uqcq_dqcq,1,h2,1,h5)
     & +(cdotpr(jz7_56_1w(:,2,1,h5),jw8_34_2w(:,2,h2))
     &  -cdotpr(jz7_56_1w(:,2,1,h5),k7561(:))
     &  *cdotpr(k7561(:),jw8_34_2w(:,2,h2))/cwmass2)/(s7561-cwmass2)
      amp(uqcq_dqcq,1,h2,1,h5)=amp(uqcq_dqcq,1,h2,1,h5)
     & +ampZ2_1728(2,h2,h5)
      enddo
      enddo
      do h5=1,2
      amp(uqcq_dqcq,1,1,1,h5)=amp(uqcq_dqcq,1,1,1,h5)
     & +ampW2_1728(2,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      amp(uqcq_dqcq,1,h2,1,h5)=amp(uqcq_dqcq,1,h2,1,h5)
     & +ampmid_1728(2,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *real(amp(uqcq_dqcq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

c      write(6,*) 'temp(2,4)',temp(2,4)
c      stop

C-----setup for (uqcq_uqsq) (2,4)->(2,3)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      amp(uqcq_uqsq,1,h2,1,h5)=
     & +cdotpr(jw8_34_2g(2,:),jz7_56_1g(:,2,h2,h5))/s7561
     & +(cdotpr(jw8_34_2z(2,:),jz7_56_1z(:,2,h2,h5))
     &  -cdotpr(jw8_34_2z(2,:),k7561(:))
     &  *cdotpr(k7561(:),jz7_56_1z(:,2,h2,h5))/czmass2)/(s7561-czmass2)
      amp(uqcq_uqsq,1,h2,1,h5)=amp(uqcq_uqsq,1,h2,1,h5)
     & +(cdotpr(jz8_56_2w(:,2,1,h5),jw7_34_1w(:,2,h2))
     &  -cdotpr(jz8_56_2w(:,2,1,h5),k7341(:))
     &  *cdotpr(k7341(:),jw7_34_1w(:,2,h2))/cwmass2)/(s7341-cwmass2)
      amp(uqcq_uqsq,1,h2,1,h5)=amp(uqcq_uqsq,1,h2,1,h5)
     & +ampZ2_2817(2,h2,h5)
      enddo
      enddo
      do h5=1,2
      amp(uqcq_uqsq,1,1,1,h5)=amp(uqcq_uqsq,1,1,1,h5)
     & +ampW2_2817(2,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      amp(uqcq_uqsq,1,h2,1,h5)=amp(uqcq_uqsq,1,h2,1,h5)
     & +ampmid_2817(2,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(2,4)=tempb(2,4)+esq**6*spinavge
     &   *real(amp(uqcq_uqsq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_uqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

C-----setup for (uquq_dquq) (2,4)->(1,4)
      ampa(uquq_dquq,:,:,:,:)=amp(uqcq_dqcq,:,:,:,:)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      ampb(uquq_dquq,1,h2,1,h5)=
     & +cdotpr(jw7_34_2g(2,:),jz8_56_1g(:,2,h2,h5))/s7342
     & +(cdotpr(jw7_34_2z(2,:),jz8_56_1z(:,2,h2,h5))
     &  -cdotpr(jw7_34_2z(2,:),k7342(:))
     &  *cdotpr(k7342(:),jz8_56_1z(:,2,h2,h5))/czmass2)/(s7342-czmass2)
      ampb(uquq_dquq,1,h2,1,h5)=ampb(uquq_dquq,1,h2,1,h5)
     & +(cdotpr(jz7_56_2w(:,2,1,h5),jw8_34_1w(:,2,h2))
     &  -cdotpr(jz7_56_2w(:,2,1,h5),k7562(:))
     &  *cdotpr(k7562(:),jw8_34_1w(:,2,h2))/cwmass2)/(s7562-cwmass2)
      ampb(uquq_dquq,1,h2,1,h5)=ampb(uquq_dquq,1,h2,1,h5)
     & +ampZ2_2718(2,h2,h5)
      enddo
      enddo
      do h5=1,2
      ampb(uquq_dquq,1,1,1,h5)=ampb(uquq_dquq,1,1,1,h5)
     & +ampW2_2718(2,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      ampb(uquq_dquq,1,h2,1,h5)=ampb(uquq_dquq,1,h2,1,h5)
     & +ampmid_2718(2,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampa(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_dquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampb(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_dquq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)-2._dp/xn*esq**6*spinavge
     &   *real(ampa(uquq_dquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_dquq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo

C-----setup for (uquq_uqdq) (2,2)->(2,1) :: not canonical order, but useful for crossings
      ampa(uquq_uqdq,:,:,:,:)=amp(uqcq_uqsq,:,:,:,:)
      if (doHO .eqv. .false.) then
      do h2=1,2
      do h5=1,2
      ampb(uquq_uqdq,1,h2,1,h5)=
     & +cdotpr(jw8_34_1g(2,:),jz7_56_2g(:,2,h2,h5))/s7562
     & +(cdotpr(jw8_34_1z(2,:),jz7_56_2z(:,2,h2,h5))
     &  -cdotpr(jw8_34_1z(2,:),k7562(:))
     &  *cdotpr(k7562(:),jz7_56_2z(:,2,h2,h5))/czmass2)/(s7562-czmass2)
      ampb(uquq_uqdq,1,h2,1,h5)=ampb(uquq_uqdq,1,h2,1,h5)
     & +(cdotpr(jz8_56_1w(:,2,1,h5),jw7_34_2w(:,2,h2))
     &  -cdotpr(jz8_56_1w(:,2,1,h5),k7342(:))
     &  *cdotpr(k7342(:),jw7_34_2w(:,2,h2))/cwmass2)/(s7342-cwmass2)
      ampb(uquq_uqdq,1,h2,1,h5)=ampb(uquq_uqdq,1,h2,1,h5)
     & +ampZ2_1827(2,h2,h5)
      enddo
      enddo
      do h5=1,2
      ampb(uquq_uqdq,1,1,1,h5)=ampb(uquq_uqdq,1,1,1,h5)
     & +ampW2_1827(2,h5)
      enddo
      endif

      do h2=1,2
      do h5=1,2
      ampb(uquq_uqdq,1,h2,1,h5)=ampb(uquq_uqdq,1,h2,1,h5)
     & +ampmid_1827(2,h2,h5)
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(2,2)=tempb(2,2)+esq**6*spinavge
     &   *real(ampa(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_uqdq,h1,h2,h3,h5)))
      tempb(2,2)=tempb(2,2)+esq**6*spinavge
     &   *real(ampb(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      tempb(2,2)=tempb(2,2)-2._dp/xn*esq**6*spinavge
     &   *real(ampa(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      endif
      enddo
      enddo
      enddo
      enddo

c--- setup for (cquq_squq) (4,2)->(3,2) :: not canonical order, but useful for crossings
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      tempb(4,2)=tempb(4,2)+esq**6*spinavge
     &   *real(ampb(uquq_uqdq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uqdq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

      if (j==1) then
      msq(2,3)=temp(2,3)
      msq(4,1)=temp(4,1)
      msq(2,1)=temp(2,1)*half
      msq(4,3)=msq(2,1)
      msq(2,4)=temp(2,4)+tempb(2,4)
      msq(2,2)=temp(2,2)
      msq(4,4)=temp(2,2)

      elseif (j==2) then
      msq(3,2)=temp(2,3)
      msq(1,4)=temp(4,1)
      msq(1,2)=temp(2,1)*half
      msq(3,4)=msq(1,2)
      msq(4,2)=temp(2,4)+tempb(2,4)

c--- qbar-qbar
      elseif (j==3) then
      msq(-1,-4)=temp(2,4)
      msq(-1,-3)=temp(2,3)
      msq(-2,-3)=tempb(2,4)
      msq(-3,-1)=temp(4,1)
      msq(-2,-1)=tempb(2,2)*half
      msq(-4,-3)=tempb(2,2)*half

c--- qbar-qbar
      elseif (j==4) then
      msq(-1,-2)=temp(2,2)*half
      msq(-3,-4)=temp(2,2)*half
      msq(-3,-2)=temp(2,4)
      msq(-1,-1)=temp(2,1)
      msq(-3,-3)=temp(2,1)
      msq(-4,-1)=tempb(2,4)
      msq(-1,-3)=msq(-1,-3)+temp(4,1)
      msq(-3,-1)=msq(-3,-1)+temp(2,3)

c--- qbar-q
      elseif (j==5) then
      msq(-3,1)=temp(2,3)
      msq(-3,3)=temp(2,1)
      msq(-1,1)=temp(2,1)
      msq(-1,3)=temp(2,3)

c--- qbar-q
      elseif (j==6) then
      msq(-4,2)=temp(2,4)
      msq(-4,4)=temp(2,2)
      msq(-3,2)=temp(2,3)+tempb(2,4)
      msq(-1,4)=msq(-3,2)
      msq(-3,4)=temp(2,1)+tempb(2,2)+temp(4,1)+tempb(4,2)
      msq(-1,2)=msq(-3,4)
      msq(-2,2)=temp(2,2)
      msq(-2,4)=temp(2,4)

c--- q-qbar
      elseif (j==7) then
      msq(1,-3)=temp(2,3)
      msq(1,-1)=temp(2,1)
      msq(2,-4)=tempb(2,4)
      msq(3,-3)=temp(2,1)
      msq(3,-1)=temp(2,3)
      msq(4,-4)=tempb(2,2)

c--- q-qbar
      elseif (j==8) then
      msq(2,-3)=tempb(2,4)+temp(2,3)
      msq(2,-2)=temp(2,2)
      msq(2,-1)=tempb(4,2)+temp(4,1)+tempb(2,2)+temp(2,1)
      msq(4,-3)=msq(2,-1)
      msq(4,-2)=temp(2,4)
      msq(4,-1)=tempb(2,4)+temp(2,3)

c      write(6,*) 'temp(2,3)',temp(2,3)
c      write(6,*) 'temp(4,1)',temp(4,1)
c      write(6,*) 'temp(2,2)',temp(2,2)
c      write(6,*) 'temp(2,1)',temp(2,1)
c      write(6,*) 'temp(2,4)',temp(2,4)
c      write(6,*) 'tempb(2,4)',tempb(2,4)
c      write(6,*) 'tempb(2,2)',tempb(2,2)
c      write(6,*) 'tempb(4,2)',tempb(4,2)
c      write(6,*)

      endif

      enddo

c--- swap elements to correct places for nwz=-1
      if (nwz .eq .-1) then
        temp(:,:)=msq(:,:)
        do j=-nf,nf
        do k=-nf,nf
        msq(j,k)=temp(-j,-k)
        enddo
        enddo
      endif

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end

