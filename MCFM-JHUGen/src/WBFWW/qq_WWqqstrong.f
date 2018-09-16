      subroutine qq_WWqqstrong(p,msq)
      implicit none
      include 'types.f'
      
c--- Author: J.M.Campbell, December 2014
c--- q(-p1)+q(-p2)->W(p3,p4)+W(p5,p6)+q(p7)+q(p8);
c--- with the t-channel exchange of a gluon.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
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
     & tempw(fn:nf,fn:nf),stat,spinavge,Colorfac,msqgg(2)
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2),ampa(nmax,2,2),ampb(nmax,2,2),cdotpr,
     & k7341(4),k8341(4),k8342(4),
     & s7341,s8341,s8342
      complex(dp)::
     & jtwo17(2,2,2,2),jtwo28(2,2,2,2),
     & jtwo18(2,2,2,2),jtwo27(2,2,2,2),
     & j7_34_1g(2,4),j8_56_2g(2,4),
     & j8_34_1g(2,4),j7_56_2g(2,4),
     & j8_34_2g(2,4),j7_56_1g(2,4)
      logical:: first
      parameter(spinavge=0.25d0,stat=0.5d0,Colorfac=V/4d0/xn**2,nfinc=4)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      data first/.true./
      save first

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
       first=.false.
       call flush(6)
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

c--- contribution from jtwoWW
      call amp2currentstrong(j1(j),j2(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,
     & jtwo17)
      call amp2currentstrong(j2(j),j1(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,
     &jtwo28)
      call amp2currentstrong(j1(j),j2(j),i3,i4,i5,i6,j8(j),j7(j),za,zb,
     & jtwo18)
      call amp2currentstrong(j2(j),j1(j),i3,i4,i5,i6,j7(j),j8(j),za,zb,
     & jtwo27)

c--- flavor-changing contributions
      call jonewstrong(j7(j),i3,i4,j1(j),za,zb,zab,j7_34_1g)
      call jonewstrong(j8(j),i5,i6,j2(j),za,zb,zab,j8_56_2g)
      call jonewstrong(j8(j),i3,i4,j1(j),za,zb,zab,j8_34_1g)
      call jonewstrong(j7(j),i5,i6,j2(j),za,zb,zab,j7_56_2g)
      call jonewstrong(j8(j),i3,i4,j2(j),za,zb,zab,j8_34_2g)
      call jonewstrong(j7(j),i5,i6,j1(j),za,zb,zab,j7_56_1g)

C-----setup for (dqcq_dqcq) 
      do h1=1,2
      do h2=1,2
      amp(dqcq_dqcq,h1,h2)=
     & +jtwo17(1,2,h1,h2)+jtwo28(2,1,h2,h1)
     
      temp(1,4)=temp(1,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(dqcq_dqcq,h1,h2)
     & *conjg(amp(dqcq_dqcq,h1,h2)))
      enddo
      enddo

C-----setup for (dqcq_uqsq) 
      amp(dqcq_uqsq,1,1)=
     & +cdotpr(j7_34_1g(1,:),j8_56_2g(2,:))/s7341

      tempw(1,4)=tempw(1,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(dqcq_uqsq,1,1)
     & *conjg(amp(dqcq_uqsq,1,1)))

C-----setup for (uqcq_uqcq) 
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(uqcq_uqcq,h1,h2)=
     & +jtwo17(2,2,h1,h2)+jtwo28(2,2,h2,h1)

      temp(2,4)=temp(2,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqcq_uqcq,h1,h2)
     & *conjg(amp(uqcq_uqcq,h1,h2)))

      enddo
      enddo


C-----setup for (dqsq_dqsq) 
      do h1=1,2
      do h2=1,2

c--- contribution from jcentre
      amp(dqsq_dqsq,h1,h2)=
     & +jtwo17(1,1,h1,h2)+jtwo28(1,1,h2,h1)

      temp(1,3)=temp(1,3)+esq**4*gsq**2*Colorfac*spinavge
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
     & +jtwo18(1,1,h1,h2)+jtwo27(1,1,h2,h1)

      temp(1,1)=temp(1,1)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(dqdq_dqdq,h1,h2)
     & *conjg(ampa(dqdq_dqdq,h1,h2)))
      temp(1,1)=temp(1,1)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(dqdq_dqdq,h1,h2)
     & *conjg(ampb(dqdq_dqdq,h1,h2)))
      if (h1 == h2) then
      temp(1,1)=temp(1,1)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
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
     & +jtwo18(2,2,h1,h2)+jtwo27(2,2,h2,h1)

      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_uquq,h1,h2)
     & *conjg(ampa(uquq_uquq,h1,h2)))
      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(uquq_uquq,h1,h2)
     & *conjg(ampb(uquq_uquq,h1,h2)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(uquq_uquq,h1,h2)
     & *conjg(ampb(uquq_uquq,h1,h2)))
      endif

      enddo
      enddo
      temp(4,4)=temp(2,2)

C-----setup for (dquq_dquq) 
c-------- ampb
      ampb(dquq_dquq,1,1)=
     & +cdotpr(j8_34_1g(1,:),j7_56_2g(2,:))/s8341

      do h1=1,2
      do h2=1,2
c-------- ampa
      ampa(dquq_dquq,h1,h2)=amp(dqcq_dqcq,h1,h2)

      temp(1,2)=temp(1,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampa(dquq_dquq,h1,h2)
     & *conjg(ampa(dquq_dquq,h1,h2)))
      temp(1,2)=temp(1,2)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(ampb(dquq_dquq,h1,h2)
     & *conjg(ampb(dquq_dquq,h1,h2)))
      if (h1 == h2) then
      temp(1,2)=temp(1,2)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
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
     & +jtwo17(2,1,h1,h2)+jtwo28(1,2,h2,h1)

      temp(2,5)=temp(2,5)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqbq_uqbq,h1,h2)
     & *conjg(amp(uqbq_uqbq,h1,h2)))

      enddo
      enddo
      temp(4,5)=temp(2,5)
      temp(2,3)=temp(2,5)

C-----setup for (uqsq_dqcq) 
      amp(uqsq_dqcq,1,1)=
     & +cdotpr(j8_34_2g(1,:),j7_56_1g(2,:))/s8342

      tempw(2,3)=tempw(2,3)+esq**4*gsq**2*Colorfac*spinavge
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

c--- 2-gluon amplitudes
      call qq2l2nuggamp(1,2,i3,i4,i5,i6,7,8,za,zb,msqgg)
      msq(1,-1)=msq(1,-1)+stat*aveqq*msqgg(1)
      msq(2,-2)=msq(2,-2)+stat*aveqq*msqgg(2)
      msq(3,-3)=msq(1,-1)+stat*aveqq*msqgg(1)
      msq(4,-4)=msq(2,-2)+stat*aveqq*msqgg(2)
      msq(5,-5)=msq(1,-1)+stat*aveqq*msqgg(1)

      call qq2l2nuggamp(2,1,i3,i4,i5,i6,7,8,za,zb,msqgg)
      msq(-1,1)=msq(-1,1)+stat*aveqq*msqgg(1)
      msq(-2,2)=msq(-2,2)+stat*aveqq*msqgg(2)
      msq(-3,3)=msq(-1,1)+stat*aveqq*msqgg(1)
      msq(-4,4)=msq(-2,2)+stat*aveqq*msqgg(2)
      msq(-5,5)=msq(-1,1)+stat*aveqq*msqgg(1)

      call qq2l2nuggamp(7,8,i3,i4,i5,i6,1,2,za,zb,msqgg)
      msq(0,0)=msq(0,0)+avegg*(3d0*msqgg(1)+2d0*msqgg(2))
c      msq(0,0)=msq(0,0)+avegg*(msqgg(1)+msqgg(2)) ! DEBUG

      call qq2l2nuggamp(7,2,i3,i4,i5,i6,1,8,za,zb,msqgg)
      msq(0,-1)=msq(0,-1)+aveqg*msqgg(1)
      msq(0,-2)=msq(0,-2)+aveqg*msqgg(2)
      msq(0,-3)=msq(0,-1)+aveqg*msqgg(1)
      msq(0,-4)=msq(0,-2)+aveqg*msqgg(2)
      msq(0,-5)=msq(0,-1)+aveqg*msqgg(1)

      call qq2l2nuggamp(7,1,i3,i4,i5,i6,2,8,za,zb,msqgg)
      msq(-1,0)=msq(-1,0)+aveqg*msqgg(1)
      msq(-2,0)=msq(-2,0)+aveqg*msqgg(2)
      msq(-3,0)=msq(-1,0)+aveqg*msqgg(1)
      msq(-4,0)=msq(-2,0)+aveqg*msqgg(2)
      msq(-5,0)=msq(-1,0)+aveqg*msqgg(1)
      
      call qq2l2nuggamp(2,8,i3,i4,i5,i6,1,7,za,zb,msqgg)
      msq(0,1)=msq(0,1)+aveqg*msqgg(1)
      msq(0,2)=msq(0,2)+aveqg*msqgg(2)
      msq(0,3)=msq(0,1)+aveqg*msqgg(1)
      msq(0,4)=msq(0,2)+aveqg*msqgg(2)
      msq(0,5)=msq(0,1)+aveqg*msqgg(1)

      call qq2l2nuggamp(1,8,i3,i4,i5,i6,2,7,za,zb,msqgg)
      msq(1,0)=msq(1,0)+aveqg*msqgg(1)
      msq(2,0)=msq(2,0)+aveqg*msqgg(2)
      msq(3,0)=msq(3,0)+aveqg*msqgg(1)
      msq(4,0)=msq(4,0)+aveqg*msqgg(2)
      msq(5,0)=msq(5,0)+aveqg*msqgg(1)

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end
      

C LocalWords:  nfinc
