      subroutine qq_WWss(p,msq)
c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W(p3,p4)+W(p5,p6)+q(p7)+q(p8);
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      include 'anom_higgs.f'
      include 'WWbits.f'
      integer:: nmax,jmax
      parameter(jmax=8,nmax=2)
      integer:: j,k,l,iq,uqcq_uqcq,uquq_uquq,
     & j1(jmax),j2(jmax),j7(jmax),j8(jmax),i3,i4,i5,i6,nfinc
      parameter(uqcq_uqcq=1,uquq_uquq=2)
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & stat,spinavge,mult
      complex(dp)::zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2),cdotpr,
     & k7341(4),k8341(4),k8342(4),k7342(4),
     & s7341,s7342,s8341,s8342
      complex(dp)::
     & jB17,jtwodiagsWpWp17,jmidWpWp17,
     & jB18,jtwodiagsWpWp18,jmidWpWp18,
     & j7_34_1z(2,4),j7_34_1g(2,4),j8_56_2z(2,4),j8_56_2g(2,4),
     & j8_34_1z(2,4),j8_34_1g(2,4),j7_56_2z(2,4),j7_56_2g(2,4),
     & j8_34_2z(2,4),j8_34_2g(2,4),j7_56_1z(2,4),j7_56_1g(2,4),
     & j8_56_1z(2,4),j8_56_1g(2,4),j7_34_2z(2,4),j7_34_2g(2,4)
      logical:: first,doHO,doBO
      parameter(spinavge=0.25_dp,stat=0.5_dp,nfinc=4)
      data j1/1,2,8,8,7,7,1,2/
      data j2/2,1,7,7,2,1,7,7/
      data j7/7,7,2,1,1,2,2,1/
      data j8/8,8,1,2,8,8,8,8/
      data first/.true./
      save j1,j2,j7,j8,first,doHO,doBO,mult
      include 'cplx.h'

      msq(:,:)=0._dp

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      if (first) then
       cwmass2=cplx2(wmass**2,-wmass*wwidth)
       czmass2=cplx2(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
c       cxw=cplx2(xw,0._dp) ! DEBUG: Madgraph comparison
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

c--- this is the MCFM ordering in process.DAT
      i3=3
      i4=4
      i5=5
      i6=6

c--- set identity of quark line based on nwz
      if (nwz == +1) then
        iq=2
      elseif (nwz == -1) then
        iq=1
      else
        write(6,*) 'Unexpected value of nwz in qq_WWss.f: ',nwz
        stop
      endif

      do j=1,jmax
      temp(:,:)=0._dp
      amp(:,:,:)=czip
      ampa(:,:,:)=czip
      ampb(:,:,:)=czip

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      k7341(:)=0.5_dp*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j7(j),:,j7(j)))
      k7342(:)=0.5_dp*(zab(j2(j),:,j2(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j7(j),:,j7(j)))
      k8341(:)=0.5_dp*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      k8342(:)=0.5_dp*(zab(j2(j),:,j2(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      s7341=cdotpr(k7341,k7341)
      s7342=cdotpr(k7342,k7342)
      s8341=cdotpr(k8341,k8341)
      s8342=cdotpr(k8342,k8342)

c--- This contribution contains Hbit and Bbit

c--- Higgs and WWWW vertex 'D' diagrams: contribution from jmidWpWp
      call ampmidWpWp
     & (j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jmidWpWp17)
      call ampmidWpWp
     & (j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jmidWpWp18)

c--- these are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then

c--- single-resonanant 'B' diagrams: contribution from ampB
      call ampBdiags(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jB17)
      call ampBdiags(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jB18)

c--- W-exchange 'C' diagrams: contribution from jtwodiagsWpWp
      call amptwodiagsWpWp
     & (j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,jtwodiagsWpWp17)
      call amptwodiagsWpWp
     & (j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,jtwodiagsWpWp18)

c--- Z/photon currents, 'A' diagrams: contribution from jtwodiagsWpWp
      call jonewWpWp(j7(j),i3,i4,j1(j),za,zb,zab,j7_34_1z,j7_34_1g)
      call jonewWpWp(j8(j),i5,i6,j2(j),za,zb,zab,j8_56_2z,j8_56_2g)
      call jonewWpWp(j8(j),i3,i4,j1(j),za,zb,zab,j8_34_1z,j8_34_1g)
      call jonewWpWp(j7(j),i5,i6,j2(j),za,zb,zab,j7_56_2z,j7_56_2g)
      call jonewWpWp(j8(j),i3,i4,j2(j),za,zb,zab,j8_34_2z,j8_34_2g)
      call jonewWpWp(j7(j),i5,i6,j1(j),za,zb,zab,j7_56_1z,j7_56_1g)
      call jonewWpWp(j8(j),i5,i6,j1(j),za,zb,zab,j8_56_1z,j8_56_1g)
      call jonewWpWp(j7(j),i3,i4,j2(j),za,zb,zab,j7_34_2z,j7_34_2g)

      else

      jB17=czip
      jB18=czip
      jtwodiagsWpWp17=czip
      jtwodiagsWpWp18=czip
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
      j8_56_1z=czip
      j8_56_1g=czip
      j7_34_2z=czip
      j7_34_2g=czip

      endif

C-----setup for (uqcq_uqcq)
      amp(uqcq_uqcq,1,1)=
     & +jB17
     & +jmidWpWp17
     & +jtwodiagsWpWp17
     & +cdotpr(j7_34_1g(iq,:),j8_56_2g(iq,:))/s7341
     & +(cdotpr(j7_34_1z(iq,:),j8_56_2z(iq,:))
     &  -cdotpr(j7_34_1z(iq,:),k7341(:))
     &  *cdotpr(k7341(:),j8_56_2z(iq,:))/czmass2)/(s7341-czmass2)
     & +cdotpr(j7_56_1g(iq,:),j8_34_2g(iq,:))/s8342
     & +(cdotpr(j7_56_1z(iq,:),j8_34_2z(iq,:))
     &  -cdotpr(j7_56_1z(iq,:),k8342(:))
     &  *cdotpr(k8342(:),j8_34_2z(iq,:))/czmass2)/(s8342-czmass2)

      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *real(amp(uqcq_uqcq,1,1)
     & *conjg(amp(uqcq_uqcq,1,1)))

C-----setup for (uquq_uquq)
c-------- ampa
      ampa(uquq_uquq,1,1)=amp(uqcq_uqcq,1,1)

c-------- ampb
      ampb(uquq_uquq,1,1)=
     & +jB18
     & +jmidWpWp18
     & +jtwodiagsWpWp18
     & +cdotpr(j8_34_1g(iq,:),j7_56_2g(iq,:))/s8341
     & +(cdotpr(j8_34_1z(iq,:),j7_56_2z(iq,:))
     &  -cdotpr(j8_34_1z(iq,:),k8341(:))
     &  *cdotpr(k8341(:),j7_56_2z(iq,:))/czmass2)/(s8341-czmass2)
     & +cdotpr(j8_56_1g(iq,:),j7_34_2g(iq,:))/s7342
     & +(cdotpr(j8_56_1z(iq,:),j7_34_2z(iq,:))
     &  -cdotpr(j8_56_1z(iq,:),k7342(:))
     &  *cdotpr(k7342(:),j7_34_2z(iq,:))/czmass2)/(s7342-czmass2)

      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampa(uquq_uquq,1,1)
     & *conjg(ampa(uquq_uquq,1,1)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     &   *real(ampb(uquq_uquq,1,1)
     & *conjg(ampb(uquq_uquq,1,1)))
      temp(2,2)=temp(2,2)-2._dp/xn*esq**6*spinavge
     &   *real(ampa(uquq_uquq,1,1)
     & *conjg(ampb(uquq_uquq,1,1)))

      temp(4,4)=temp(2,2)

c--- fill matrix elements
      if (j == 1) then

      do k=1,nfinc
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nfinc
      msq(k,l)=temp(k,l)
      enddo
      enddo

      elseif (j==2) then
      do k=1,nfinc
      do l=k+1,nfinc
      msq(l,k)=temp(k,l)
      enddo
      enddo

      elseif (j==3) then
      msq(-3,-1)=temp(2,4)
      msq(-1,-1)=temp(2,2)*stat
      msq(-3,-3)=temp(4,4)*stat

      elseif (j==4) then
      msq(-1,-3)=temp(2,4)

c--- qbar-q
      elseif (j==5) then
      msq(-1,4)=temp(2,4)
      msq(-3,2)=temp(2,4)
      msq(-1,2)=temp(2,2) ! d~u -> u~d

c--- q-qbar
      elseif (j==6) then
      msq(2,-1)=temp(2,2) ! ud~ -> u~d
      msq(4,-1)=temp(2,4)
      msq(2,-3)=temp(2,4)

c--- q-qbar extra pieces
      elseif (j==7) then
      msq(2,-1)=msq(2,-1)+temp(2,4) ! ud~ -> c~s
      msq(4,-3)=msq(2,-1)

c--- qbar-q extra pieces
      elseif (j==8) then
      msq(-1,2)=msq(-1,2)+temp(2,4) ! d~u -> c~s
      msq(-3,4)=msq(-1,2)

      endif

      enddo

c--- above assignments assume W+W-; re-assign for W-W-
      if (nwz == -1) then
        temp(:,:)=msq(:,:)
        msq(:,:)=0._dp
        msq(-4,-4)=temp(-3,-3)
        msq(-4,-2)=temp(-3,-1)
        msq(-4,1)=temp(-3,2)
        msq(-4,3)=temp(-3,4)
        msq(-2,-4)=temp(-1,-3)
        msq(-2,-2)=temp(-1,-1)
        msq(-2,1)=temp(-1,2)
        msq(-2,3)=temp(-1,4)
        msq(1,-4)=temp(2,-3)
        msq(1,-2)=temp(2,-1)
        msq(1,1)=temp(2,2)
        msq(1,3)=temp(2,4)
        msq(3,-4)=temp(4,-3)
        msq(3,-2)=temp(4,-1)
        msq(3,1)=temp(4,2)
        msq(3,3)=temp(4,4)
      endif

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end

