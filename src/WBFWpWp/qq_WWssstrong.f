      subroutine qq_WWssstrong(p,msq)
      implicit none
c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W+(p3,p4)+W+(p5,p6)+q(p7)+q(p8);
c--- with the t-channel exchange of a gluon.
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      integer nmax,jmax
      parameter(jmax=8,nmax=2)
      integer j,k,l,iq,uqcq_uqcq,uquq_uquq,
     & j1(jmax),j2(jmax),j7(jmax),j8(jmax),i3,i4,i5,i6,nfinc
      parameter(uqcq_uqcq=1,uquq_uquq=2)
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & stat,spinavge,Colorfac
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),cdotpr,
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2),
     & k7341(4),k8341(4),k8342(4),k7342(4),
     & s7341,s7342,s8341,s8342
      double complex
     & j7_34_1g(2,4),j8_56_2g(2,4),j8_34_1g(2,4),j7_56_2g(2,4),
     & j8_34_2g(2,4),j7_56_1g(2,4),j8_56_1g(2,4),j7_34_2g(2,4)
      logical first
      parameter(spinavge=0.25d0,stat=0.5d0,Colorfac=V/4d0/xn**2,nfinc=4)
      data j1/1,2,8,8,7,7,1,2/
      data j2/2,1,7,7,2,1,7,7/
      data j7/7,7,2,1,1,2,2,1/
      data j8/8,8,1,2,8,8,8,8/
      data first/.true./
      save j1,j2,j7,j8,first

      msq(:,:)=0d0

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      cwmass2=dcmplx(wmass**2,0d0)
      czmass2=dcmplx(zmass**2,0d0)
      cxw=dcmplx(xw,0d0)
      
c--- this is the MCFM ordering in process.DAT
      i3=3
      i4=4
      i5=5
      i6=6

c--- set identity of quark line based on nwz
      if (nwz .eq. +1) then
        iq=2
      elseif (nwz .eq. -1) then
        iq=1
      else
        write(6,*) 'Unexpected value of nwz in qq_WWss.f: ',nwz
        stop
      endif

      do j=1,jmax
      temp(:,:)=0d0
      amp(:,:,:)=czip
      ampa(:,:,:)=czip
      ampb(:,:,:)=czip

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      k7341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j7(j),:,j7(j)))
      k7342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j7(j),:,j7(j)))
      k8341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      k8342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(i3,:,i3)
     & +zab(i4,:,i4)+zab(j8(j),:,j8(j)))
      s7341=cdotpr(k7341,k7341)
      s7342=cdotpr(k7342,k7342)
      s8341=cdotpr(k8341,k8341)
      s8342=cdotpr(k8342,k8342)

c--- Z/photon currents, 'A' diagrams: contribution from jtwodiagsWpWp
      call jonewstrong(j7(j),i3,i4,j1(j),za,zb,zab,j7_34_1g)
      call jonewstrong(j8(j),i5,i6,j2(j),za,zb,zab,j8_56_2g)
      call jonewstrong(j8(j),i3,i4,j1(j),za,zb,zab,j8_34_1g)
      call jonewstrong(j7(j),i5,i6,j2(j),za,zb,zab,j7_56_2g)
      call jonewstrong(j8(j),i3,i4,j2(j),za,zb,zab,j8_34_2g)
      call jonewstrong(j7(j),i5,i6,j1(j),za,zb,zab,j7_56_1g)
      call jonewstrong(j8(j),i5,i6,j1(j),za,zb,zab,j8_56_1g)
      call jonewstrong(j7(j),i3,i4,j2(j),za,zb,zab,j7_34_2g)
      
C-----setup for (uqcq_uqcq) 
      amp(uqcq_uqcq,1,1)=
     & +cdotpr(j7_34_1g(iq,:),j8_56_2g(iq,:))/s7341
     & +cdotpr(j7_56_1g(iq,:),j8_34_2g(iq,:))/s8342
      temp(2,4)=temp(2,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *dble(amp(uqcq_uqcq,1,1)
     & *dconjg(amp(uqcq_uqcq,1,1)))
      
C-----setup for (uquq_uquq) 
c-------- ampa
      ampa(uquq_uquq,1,1)=amp(uqcq_uqcq,1,1)

c-------- ampb
      ampb(uquq_uquq,1,1)=
     & +cdotpr(j8_34_1g(iq,:),j7_56_2g(iq,:))/s8341
     & +cdotpr(j8_56_1g(iq,:),j7_34_2g(iq,:))/s7342

      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *dble(ampa(uquq_uquq,1,1)
     & *dconjg(ampa(uquq_uquq,1,1)))
      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *dble(ampb(uquq_uquq,1,1)
     & *dconjg(ampb(uquq_uquq,1,1)))
      temp(2,2)=temp(2,2)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
     &   *dble(ampa(uquq_uquq,1,1)
     & *dconjg(ampb(uquq_uquq,1,1)))

      temp(4,4)=temp(2,2)

c--- fill matrix elements
      if (j .eq. 1) then

      do k=1,nfinc
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nfinc
      msq(k,l)=temp(k,l)
      enddo
      enddo

      elseif (j.eq.2) then
      do k=1,nfinc
      do l=k+1,nfinc
      msq(l,k)=temp(k,l)
      enddo
      enddo

      elseif (j.eq.3) then
      msq(-3,-1)=temp(2,4)
      msq(-1,-1)=temp(2,2)*stat
      msq(-3,-3)=temp(4,4)*stat

      elseif (j.eq.4) then
      msq(-1,-3)=temp(2,4)

c--- qbar-q
      elseif (j.eq.5) then
      msq(-1,4)=temp(2,4)
      msq(-3,2)=temp(2,4)
      msq(-1,2)=temp(2,2) ! d~u -> u~d
      
c--- q-qbar
      elseif (j.eq.6) then
      msq(2,-1)=temp(2,2) ! ud~ -> u~d
      msq(4,-1)=temp(2,4)
      msq(2,-3)=temp(2,4)

c--- q-qbar extra pieces
      elseif (j.eq.7) then
      msq(2,-1)=msq(2,-1)+temp(2,4) ! ud~ -> c~s
      msq(4,-3)=msq(2,-1)

c--- qbar-q extra pieces
      elseif (j.eq.8) then
      msq(-1,2)=msq(-1,2)+temp(2,4) ! d~u -> c~s
      msq(-3,4)=msq(-1,2)

      endif
      
      enddo

c--- above assignments assume W+W-; re-assign for W-W-
      if (nwz .eq. -1) then
        temp(:,:)=msq(:,:)
        msq(:,:)=0d0
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
      
