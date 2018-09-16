      subroutine qq_ZZqqstrong(p,msq)
      implicit none
      include 'types.f'
      
c--- Author: R.K. Ellis, October 2014
c--- q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
c--- with the t-channel exchange of a gluon.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'first.f'
      integer:: nmax,jmax
      parameter(jmax=12,nmax=10)
      integer:: j,k,l,i1,i2,i3,i4,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq
c     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6)
c     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer:: h1,h2,h3,h5
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & stat,spinavge,Colorfac,t4,
     & s17,s28,s18,s27,s7341,s7561,s7342,s7562,
     & msqgg(2)
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j7_1(4,2),j7_2(4,2),j8_1(4,2),j8_2(4,2),cdotpr,
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
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2)
      parameter(spinavge=0.25d0,stat=0.5d0,Colorfac=V/4d0/xn**2)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
c--- Begin statement functions
      t4(i1,i2,i3,i4)=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4) 
     & +s(i2,i3)+s(i2,i4)+s(i3,i4) 
c--- End statement functions

      msq(:,:)=0d0

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)
      if (first) then
       cwmass2=cplx2(wmass**2,-wmass*wwidth)
       czmass2=cplx2(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
c       cxw=cplx2(xw,0d0) ! DEBUG: Madgraph comparison
c       gsq=(-1.2177157847767197d0)**2 ! DEBUG: Madgraph comparison
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
      endif
      
C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      do j=1,jmax

      temp(:,:)=0d0
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

      s17=s(j1(j),j7(j))
      s28=s(j2(j),j8(j))
      s27=s(j2(j),j7(j))
      s18=s(j1(j),j8(j))
      s7341=t4(j7(j),3,4,j1(j))
      s7342=t4(j7(j),3,4,j2(j))
      s7561=t4(j7(j),5,6,j1(j))
      s7562=t4(j7(j),5,6,j2(j))

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

      amp=czip

C-----setup for (uqbq_uqbq) (2,5)->(2,5)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C---one-one currents
      amp(uqbq_uqbq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))/s7561

C---two-one currents
      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,1,h2,h3,h5))/s17

      temp(2,5)=temp(2,5)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *conjg(amp(uqbq_uqbq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(2,3)=temp(2,5)
      temp(4,5)=temp(2,5)
C--------------------------------------------------------------------------
C-----setup for (uqcq_uqcq) (2,4)->(2,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

C---one-one currents
      amp(uqcq_uqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))/s7561
      
C---two-one currents
      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))/s17

      temp(2,4)=temp(2,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *conjg(amp(uqcq_uqcq,h1,h2,h3,h5)))

      enddo
      enddo
      enddo
      enddo


C-----setup for (dqcq_dqcq) (1,4)-->(1,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(dqcq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))/s7561

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))/s17

      temp(1,4)=temp(1,4)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *conjg(amp(dqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(1,2)=temp(1,4)
      temp(3,4)=temp(1,4)

C-----setup for (dqsq_dqsq) (1,3)-->(1,3)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      amp(dqsq_dqsq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))/s7341
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))/s7561

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))/s28
     & +cdotpr(j8_3456_2(:,1,h2,h3,h5),j7_1(:,h1))/s17

      temp(1,3)=temp(1,3)+esq**4*gsq**2*Colorfac*spinavge
     &   *real(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *conjg(amp(dqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

      if ((j==2).or.(j==4).or.(j==6).or.(j==8)
     & .or.(j==10).or.(j==12)) go to 100
C-----setup for ((uquq_uquq)  (2,2)-->(2,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C-----------------ampa
      ampa(uquq_uquq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)

C-----------------ampb
      ampb(uquq_uquq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))/s7562
     & +cdotpr(j8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))/s7342

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,2,h1,h3,h5),j7_2(:,h2))/s27
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,2,h2,h3,h5))/s18

      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     & *real(ampa(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampa(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**4*gsq**2*Colorfac*spinavge
     & *real(ampb(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(2,2)=temp(2,2)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
     & *real(ampa(uquq_uquq,h1,h2,h3,h5)
     & *conjg(ampb(uquq_uquq,h1,h2,h3,h5)))
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
      ampb(dqdq_dqdq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))/s7562
     & +cdotpr(j8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))/s7342

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,1,h1,h3,h5),j7_2(:,h2))/s27
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,1,h2,h3,h5))/s18

      temp(1,1)=temp(1,1)+esq**4*gsq**2*Colorfac*spinavge
     & *real(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampa(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+esq**4*gsq**2*Colorfac*spinavge
     & *real(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 == h2) then
      temp(1,1)=temp(1,1)+2d0/xn*esq**4*gsq**2*Colorfac*spinavge
     & *real(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *conjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)
 100  continue

      if (j==1) then
      do k=1,nf
      msq(k,k)=temp(k,k)*stat
      do l=k+1,nf
      msq(k,l)=temp(k,l)
      enddo
      enddo

      elseif (j==2) then
      do k=1,nf
      do l=k+1,nf
      msq(l,k)=temp(k,l)
      enddo
      enddo

      elseif (j==3) then
      do k=-nf,-1
      msq(k,k)=temp(-k,-k)*stat
      do l=k+1,-1
      msq(k,l)=temp(-l,-k)
      enddo
      enddo

      elseif (j==4) then
      do k=-nf,-1
      do l=k+1,-1
      msq(l,k)=temp(-l,-k)
      enddo
      enddo

c--- qbar-q
      elseif (j==5) then
      do k=-nf,-1
      msq(k,-k)=temp(-k,-k)
      do l=1,nf
      if (abs(k) < abs(l)) then
      msq(k,l)=temp(-k,l)
      endif
      enddo
      enddo
      
c--- qbar-q
      elseif (j==6) then
      do k=-nf,-1
      do l=1,nf
      if (abs(k) > abs(l)) then
      msq(k,l)=temp(l,-k)
      endif
      enddo
      enddo

c--- q-qbar
      elseif (j==7) then
      do k=-nf,-1
      msq(-k,k)=temp(-k,-k)
      do l=1,nf
      if (abs(k) < abs(l)) then
      msq(l,k)=temp(-k,l)
      endif
      enddo
      enddo

c--- q-qbar
      elseif (j==8) then
      do k=-nf,-1
      do l=-nf,-1
      if (abs(k) < abs(l)) then
      msq(-k,l)=temp(-k,-l)
      endif
      enddo
      enddo
      
c--- q-qbar extra pieces
      elseif (j==9) then
      do k=1,nf
      do l=1,nf
      if (k < l) then
      msq(k,-k)=msq(k,-k)+temp(k,l)
      endif
      enddo
      enddo
      msq(3,-4)=msq(1,-2)
      msq(4,-3)=msq(2,-1)
c--- q-qbar extra pieces
      elseif (j==10) then
      do k=1,nf
      do l=1,nf
      if (k > l) then
      msq(k,-k)=msq(k,-k)+temp(l,k)
      endif
      enddo
      enddo
 
c--- qbar-q extra pieces
      elseif (j==11) then
      do k=1,nf
      do l=1,nf
      if (k < l) then
      msq(-k,k)=msq(-k,k)+temp(k,l)
      endif
      enddo
      enddo
      msq(-4,3)=msq(-2,1)
      msq(-3,4)=msq(-1,2)

c--- qbar-q extra pieces
      elseif (j==12) then
      do k=1,nf
      do l=1,nf
      if (k > l) then
      msq(-k,k)=msq(-k,k)+temp(l,k)
      endif
      enddo
      enddo
  
      endif

      enddo

c--- 2-gluon amplitudes
      call qq4lggampf(1,2,3,4,5,6,7,8,3,4,za,zb,msqgg)
      msq(1,-1)=msq(1,-1)+stat*aveqq*msqgg(1)
      msq(2,-2)=msq(2,-2)+stat*aveqq*msqgg(2)
      msq(3,-3)=msq(1,-1)+stat*aveqq*msqgg(1)
      msq(4,-4)=msq(2,-2)+stat*aveqq*msqgg(2)
      msq(5,-5)=msq(1,-1)+stat*aveqq*msqgg(1)

      call qq4lggampf(2,1,3,4,5,6,7,8,3,4,za,zb,msqgg)
      msq(-1,1)=msq(-1,1)+stat*aveqq*msqgg(1)
      msq(-2,2)=msq(-2,2)+stat*aveqq*msqgg(2)
      msq(-3,3)=msq(-1,1)+stat*aveqq*msqgg(1)
      msq(-4,4)=msq(-2,2)+stat*aveqq*msqgg(2)
      msq(-5,5)=msq(-1,1)+stat*aveqq*msqgg(1)

      call qq4lggampf(7,8,3,4,5,6,1,2,3,4,za,zb,msqgg)
      msq(0,0)=msq(0,0)+avegg*(3d0*msqgg(1)+2d0*msqgg(2))

      call qq4lggampf(7,2,3,4,5,6,1,8,3,4,za,zb,msqgg)
      msq(0,-1)=msq(0,-1)+aveqg*msqgg(1)
      msq(0,-2)=msq(0,-2)+aveqg*msqgg(2)
      msq(0,-3)=msq(0,-1)+aveqg*msqgg(1)
      msq(0,-4)=msq(0,-2)+aveqg*msqgg(2)
      msq(0,-5)=msq(0,-1)+aveqg*msqgg(1)

      call qq4lggampf(7,1,3,4,5,6,2,8,3,4,za,zb,msqgg)
      msq(-1,0)=msq(-1,0)+aveqg*msqgg(1)
      msq(-2,0)=msq(-2,0)+aveqg*msqgg(2)
      msq(-3,0)=msq(-1,0)+aveqg*msqgg(1)
      msq(-4,0)=msq(-2,0)+aveqg*msqgg(2)
      msq(-5,0)=msq(-1,0)+aveqg*msqgg(1)
      
      call qq4lggampf(2,8,3,4,5,6,1,7,3,4,za,zb,msqgg)
      msq(0,1)=msq(0,1)+aveqg*msqgg(1)
      msq(0,2)=msq(0,2)+aveqg*msqgg(2)
      msq(0,3)=msq(0,1)+aveqg*msqgg(1)
      msq(0,4)=msq(0,2)+aveqg*msqgg(2)
      msq(0,5)=msq(0,1)+aveqg*msqgg(1)

      call qq4lggampf(1,8,3,4,5,6,2,7,3,4,za,zb,msqgg)
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
      
